#include "DNNInferencer_blas.H"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <yaml-cpp/yaml.h>

DNNInferencer_blas::DNNInferencer_blas() {
    char* env_tmp = getenv("DNN_BATCH_SIZE");
    if(env_tmp == NULL){
        this->batch_size_ = 16384;
    }else{
        this->batch_size_ = std::atol(env_tmp);;
    }
}                

DNNInferencer_blas::~DNNInferencer_blas() {
    for(int i = 0;i < model0_.size(); ++i){
        delete model0_[i];
    }
    model0_.clear();
    for(int i = 0;i < model1_.size(); ++i){
        delete model1_[i];
    }
    model1_.clear();
    for(int i = 0;i < model2_.size(); ++i){
        delete model2_[i];
    }
    model2_.clear();
}


void DNNInferencer_blas::load_models(const std::string dir){
    // mpi 
    int flag_mpi_init;
    MPI_Initialized(&flag_mpi_init);
    if(!flag_mpi_init){
        std::cerr << "DNNInferencer_blas::load_models : MPI is not initialized" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int mpirank;
    int mpisize;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

    int32_t count;
    char* buffer;
    std::string setting0_str;

    if(mpirank == 0){
        std::ifstream fin(dir + "/0/setting.yaml");
        if (!fin) {
            std::cerr << "open setting file error , setting path : " << dir + "/0/setting.yaml" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        std::ostringstream oss;
        oss << fin.rdbuf();
        fin.close();
        setting0_str = oss.str();
        count = setting0_str.size();
        buffer = new char[count];
        std::copy(setting0_str.begin(), setting0_str.end(), buffer);
    }

    MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(mpirank != 0){
        buffer = new char[count];
    }

    MPI_Bcast(buffer, count, MPI_CHAR, 0, MPI_COMM_WORLD);

    if(mpirank != 0){
        setting0_str = std::string(buffer, count);
    }

    delete[] buffer;

    // init model
    YAML::Node setting0 = YAML::Load(setting0_str);
    YAML::Node layersNode0 = setting0["layers"];
    for(size_t i = 0; i < layersNode0.size(); ++i){
        layers_.push_back(layersNode0[i].as<int64_t>());
    }

    YAML::Node modelNode0 = setting0["model"];
    for(size_t i = 0; i < modelNode0.size(); ++i){
        std::string layerType = modelNode0[i]["layer"]["type"].as<std::string>();
        // std::string weight_path = modelNode0[i]["layer"]["weight_path"].as<std::string>();
        // std::string bias_path = modelNode0[i]["layer"]["bias_path"].as<std::string>();
        int64_t in_features = modelNode0[i]["layer"]["in_features"].as<int64_t>();
        int64_t out_features = modelNode0[i]["layer"]["out_features"].as<int64_t>();
        if(layerType == "LinearGELU"){
            model0_.push_back(new LinearGELU<float>(in_features, out_features));
        }else if(layerType == "Linear"){
            model0_.push_back(new Linear<float>(in_features, out_features));
        }else{
            assert(false);
        }
    }

    YAML::Node modelNode1 = setting0["model"];
    for(size_t i = 0; i < modelNode1.size(); ++i){
        std::string layerType = modelNode1[i]["layer"]["type"].as<std::string>();
        // std::string weight_path = modelNode1[i]["layer"]["weight_path"].as<std::string>();
        // std::string bias_path = modelNode1[i]["layer"]["bias_path"].as<std::string>();
        int64_t in_features = modelNode1[i]["layer"]["in_features"].as<int64_t>();
        int64_t out_features = modelNode1[i]["layer"]["out_features"].as<int64_t>();
        if(layerType == "LinearGELU"){
            model1_.push_back(new LinearGELU<float>(in_features, out_features));
        }else if(layerType == "Linear"){
            model1_.push_back(new Linear<float>(in_features, out_features));
        }else{
            assert(false);
        }
    }
    
    YAML::Node modelNode2 = setting0["model"];
    for(size_t i = 0; i < modelNode2.size(); ++i){
        std::string layerType = modelNode2[i]["layer"]["type"].as<std::string>();
        // std::string weight_path = modelNode2[i]["layer"]["weight_path"].as<std::string>();
        // std::string bias_path = modelNode2[i]["layer"]["bias_path"].as<std::string>();
        int64_t in_features = modelNode2[i]["layer"]["in_features"].as<int64_t>();
        int64_t out_features = modelNode2[i]["layer"]["out_features"].as<int64_t>();
        if(layerType == "LinearGELU"){
            model2_.push_back(new LinearGELU<float>(in_features, out_features));
        }else if(layerType == "Linear"){
            model2_.push_back(new Linear<float>(in_features, out_features));
        }else{
            assert(false);
        }
    }

    // load parameters
    for(int i = 0; i < model0_.size(); ++i){
        model0_[i]->load_parameters(dir+"/0", i);
    }
    for(int i = 0; i < model1_.size(); ++i){
        model1_[i]->load_parameters(dir+"/1", i);
    }
    for(int i = 0; i < model2_.size(); ++i){
        model2_[i]->load_parameters(dir+"/2", i);
    }

    buffer_alloced_ = true;
    FLOPs_per_sample_ = 0;
    for(size_t i = 1; i < layers_.size(); ++i){
        output_buffer_.emplace_back(std::vector<float>(batch_size_ * layers_[i]));
        FLOPs_per_sample_ += 2.0 * layers_[i - 1] * layers_[i];
    }
}

void DNNInferencer_blas::Inference_multiDNNs(
    const std::vector<float>& input0, std::vector<double>& output0, int64_t input_count0,
    const std::vector<float>& input1, std::vector<double>& output1, int64_t input_count1,
    const std::vector<float>& input2, std::vector<double>& output2, int64_t input_count2
){
    double dnn_infer_start = MPI_Wtime();

    if(input_count0 > 0){
        output0.resize(input_count0 * output_dim());

        for(int64_t sample_start = 0; sample_start < input_count0; sample_start += batch_size_){
            int64_t sample_end = std::min(input_count0, sample_start + batch_size_);
            int64_t sample_len = sample_end - sample_start;
            std::vector<Tensor<float>> tensor_list;
            tensor_list.emplace_back(Tensor<float>({sample_len, layers_[0]}, const_cast<float*>(input0.data()) + sample_start * input_dim()));
            for(size_t i = 1; i < layers_.size(); ++i){
                tensor_list.emplace_back(Tensor<float>({sample_len, layers_[i]}, output_buffer_[i - 1].data()));
            }
            for(size_t i = 0; i < model0_.size(); ++i){
                model0_[i]->forward(tensor_list[i], tensor_list[i+1]);
            }
            Tensor<float>& last_tensor = tensor_list.back();

            double* __restrict__ output0_ptr = output0.data() + sample_start * output_dim();
            const float* const __restrict__ last_tensor_ptr = last_tensor.data();

            for(int i = 0; i < last_tensor.element_num(); ++i){
                output0_ptr[i] = last_tensor_ptr[i];
            }
        }
    }

    if(input_count1 > 0){
        output1.resize(input_count1 * output_dim());

        for(int64_t sample_start = 0; sample_start < input_count1; sample_start += batch_size_){
            int64_t sample_end = std::min(input_count1, sample_start + batch_size_);
            int64_t sample_len = sample_end - sample_start;

            std::vector<Tensor<float>> tensor_list;
            tensor_list.emplace_back(Tensor<float>({sample_len, layers_[0]}, const_cast<float*>(input1.data()) + sample_start * input_dim()));
            for(size_t i = 1; i < layers_.size(); ++i){
                tensor_list.emplace_back(Tensor<float>({sample_len, layers_[i]}, output_buffer_[i - 1].data()));
            }

            for(size_t i = 0; i < model1_.size(); ++i){
                model1_[i]->forward(tensor_list[i], tensor_list[i+1]);
            }

            Tensor<float>& last_tensor = tensor_list.back();

            double* __restrict__ output1_ptr = output1.data() + sample_start * output_dim();
            const float* const __restrict__ last_tensor_ptr = last_tensor.data();

            for(int i = 0; i < last_tensor.element_num(); ++i){
                output1_ptr[i] = last_tensor_ptr[i];
            }
        }
    }

    if(input_count2 > 0){
        output2.resize(input_count2 * output_dim());

        for(int64_t sample_start = 0; sample_start < input_count2; sample_start += batch_size_){
            int64_t sample_end = std::min(input_count2, sample_start + batch_size_);
            int64_t sample_len = sample_end - sample_start;

            std::vector<Tensor<float>> tensor_list;
            tensor_list.emplace_back(Tensor<float>({sample_len, layers_[0]}, const_cast<float*>(input2.data()) + sample_start * input_dim()));
            for(size_t i = 1; i < layers_.size(); ++i){
                tensor_list.emplace_back(Tensor<float>({sample_len, layers_[i]}, output_buffer_[i - 1].data()));
            }

            for(size_t i = 0; i < model2_.size(); ++i){
                model2_[i]->forward(tensor_list[i], tensor_list[i+1]);
            }

            Tensor<float>& last_tensor = tensor_list.back();

            double* __restrict__ output2_ptr = output2.data() + sample_start * output_dim();
            const float* const __restrict__ last_tensor_ptr = last_tensor.data();

            for(int i = 0; i < last_tensor.element_num(); ++i){
                output2_ptr[i] = last_tensor_ptr[i];
            }
        }
    }

    double dnn_infer_end = MPI_Wtime();
    double dnn_infer_time = dnn_infer_end - dnn_infer_start;
    double FLOPs = (input_count0 + input_count1 + input_count2) * FLOPs_per_sample_;
    
#ifdef _OPENMP
    int num_threads = omp_get_max_threads();
#else
    int num_threads = 1;
#endif
    double theoretical_peak = 3.3792 / 48. * 2. * num_threads;
    double FLOPS = FLOPs / dnn_infer_time;
    double TFLOPS = FLOPS * 1e-12;
    double peak = TFLOPS * 100. / theoretical_peak;

    int mpirank;
    int flag_mpi_init;
    MPI_Initialized(&flag_mpi_init);

    if(flag_mpi_init) MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

    if(mpirank == 0){
        std::cout << "Inference Performance ---------------" << std::endl;
        std::cout << "samples : " << (input_count0 + input_count1 + input_count2) << std::endl;
        std::cout << "batch size : " << batch_size_ << std::endl;
        std::cout << "Time : " << dnn_infer_time << std::endl;
        std::cout << "FLOPS : " << FLOPs << std::endl;
        std::cout << "TFLOPS : " << TFLOPS << std::endl;
        std::cout << "Peak : " << peak << std::endl;
        std::cout << "-------------------------------------" << std::endl;
    } 
}



