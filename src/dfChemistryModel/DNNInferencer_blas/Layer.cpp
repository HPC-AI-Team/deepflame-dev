#include "Layer.H"
#include <cmath>
#include <sstream>
#include <fstream>
#include <cassert>
#include <iostream>
#include <mpi.h>
    
#ifdef __cplusplus
extern "C" {
#endif

extern void sgemm_(char *transa, char *transb, int *m, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);
extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

#ifdef __cplusplus
}
#endif

inline void gemm(char *transa, char *transb, int *m, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc){
    sgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void gemm(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc){
    dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

void read_float_data(float* data, int64_t len, std::string filepath){
    std::ifstream fin(filepath, std::ios::binary);
    if(!fin.is_open()){
        std::cerr << "open file error : " << filepath << std::endl << std::flush;
        abort();
    }
    fin.read(reinterpret_cast<char*>(data), len * sizeof(float));
    fin.close();
}

template<>
void Linear<float>::load_parameters(const std::string& dir, int64_t layer_id){
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

    if(mpirank == 0){
        std::stringstream ss1,ss2;
        ss1 << dir << "/" << "linear_" << layer_id << "_weights_rowmajor_" << in_features_ << "_" << out_features_ << ".data";
        std::string weights_path = ss1.str();
        ss2 << dir << "/" << "linear_" << layer_id << "_bias_" << out_features_ << ".data";
        std::string bias_path = ss2.str();

        read_float_data(weights_.data(), weights_.element_num(), weights_path);
        read_float_data(bias_.data(), bias_.element_num(), bias_path);
    }

    MPI_Bcast(weights_.data(), weights_.element_num(), MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(bias_.data(), bias_.element_num(), MPI_FLOAT, 0, MPI_COMM_WORLD);
}

template<>
void LinearGELU<float>::load_parameters(const std::string& dir, int64_t layer_id){
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

    if(mpirank == 0){
        std::stringstream ss1,ss2;
        ss1 << dir << "/" << "linear_" << layer_id << "_weights_rowmajor_" << in_features_ << "_" << out_features_ << ".data";
        std::string weights_path = ss1.str();
        ss2 << dir << "/" << "linear_" << layer_id << "_bias_" << out_features_ << ".data";
        std::string bias_path = ss2.str();

        read_float_data(weights_.data(), weights_.element_num(), weights_path);
        read_float_data(bias_.data(), bias_.element_num(), bias_path);
    }

    MPI_Bcast(weights_.data(), weights_.element_num(), MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(bias_.data(), bias_.element_num(), MPI_FLOAT, 0, MPI_COMM_WORLD);
}

template<>
void Linear<double>::load_parameters(const std::string& dir, int64_t layer_id){
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

    if(mpirank == 0){
        std::stringstream ss1,ss2;
        ss1 << dir << "/" << "linear_" << layer_id << "_weights_rowmajor_" << in_features_ << "_" << out_features_ << ".data";
        std::string weights_path = ss1.str();
        ss2 << dir << "/" << "linear_" << layer_id << "_bias_" << out_features_ << ".data";
        std::string bias_path = ss2.str();

        float* weight_tmp = new float[weights_.element_num()];
        float* bias_tmp = new float[bias_.element_num()];

        read_float_data(weight_tmp, weights_.element_num(), weights_path);
        read_float_data(bias_tmp, bias_.element_num(), bias_path);

        double* weights_data = weights_.data();
        double* bias_data = bias_.data();
        for(int i = 0; i < weights_.element_num(); ++i)
            weights_data[i] = static_cast<double>(weight_tmp[i]);
        for(int i = 0; i < bias_.element_num(); ++i)
            bias_data[i] = static_cast<double>(bias_tmp[i]);

        delete [] weight_tmp;
        delete [] bias_tmp;
    }

    MPI_Bcast(weights_.data(), weights_.element_num(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(bias_.data(), bias_.element_num(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

template<>
void LinearGELU<double>::load_parameters(const std::string& dir, int64_t layer_id){
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

    if(mpirank == 0){
        std::stringstream ss1,ss2;
        ss1 << dir << "/" << "linear_" << layer_id << "_weights_rowmajor_" << in_features_ << "_" << out_features_ << ".data";
        std::string weights_path = ss1.str();
        ss2 << dir << "/" << "linear_" << layer_id << "_bias_" << out_features_ << ".data";
        std::string bias_path = ss2.str();
        
        float* weight_tmp = new float[weights_.element_num()];
        float* bias_tmp = new float[bias_.element_num()];

        read_float_data(weight_tmp, weights_.element_num(), weights_path);
        read_float_data(bias_tmp, bias_.element_num(), bias_path);

        double* weights_data = weights_.data();
        double* bias_data = bias_.data();
        for(int i = 0; i < weights_.element_num(); ++i)
            weights_data[i] = static_cast<double>(weight_tmp[i]);
        for(int i = 0; i < bias_.element_num(); ++i)
            bias_data[i] = static_cast<double>(bias_tmp[i]);

        delete [] weight_tmp;
        delete [] bias_tmp;
    }

    MPI_Bcast(weights_.data(), weights_.element_num(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(bias_.data(), bias_.element_num(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


void gelu_navie(int64_t len, float* data){
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int64_t i = 0; i < len; ++i){
        float x = data[i];
        data[i] = 0.5f * x * (1.f + tanhf(sqrtf(2.f / M_PI) * (x + 0.044715f * powf(x, 3.f))));
    }
}

void gelu_navie(int64_t len, double* data){
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int64_t i = 0; i < len; ++i){
        double x = data[i];
        data[i] = 0.5 * x * (1. + tanh(sqrt(2. / M_PI) * (x + 0.044715 * pow(x, 3.))));
    }
}

// inline float tanh_exp(float x){
//     return 1.f - 2.f / (expf(2.f * x) + 1.f);
// }

// inline double tanh_exp(double x){
//     return 1. - 2. / (exp(2. * x) + 1.);
// }

inline float tanh_exp(float x){
    if(x > 8.f) return 1.f;
    if(x < -8.f) return -1.f;
    return 1.f - 2.f / (expf(2.f * x) + 1.f);
}

inline double tanh_exp(double x){
    if(x > 8.) return 1.;
    if(x < -8.) return -1.;
    return 1. - 2. / (exp(2. * x) + 1.);
}

void gelu_exp(int64_t len, float* data){
    const float const_1 = sqrtf(2.f / M_PI);
    const float const_2 = 0.044715f;
    const float one = 1.f;
    const float half = 0.5;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int64_t i = 0; i < len; ++i){
        float x = data[i];
        data[i] = half * x * (one + tanh_exp(const_1 * (x + const_2 * x * x * x)));
    }
}

void gelu_exp(int64_t len, double* data){
    const double const_1 = sqrt(2. / M_PI);
    const double const_2 = 0.044715;
    const double one = 1.;
    const double half = 0.5;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int64_t i = 0; i < len; ++i){
        double x = data[i];
        data[i] = half * x * (one + tanh_exp(const_1 * (x + const_2 * x * x * x)));
    }
}

void bias_naive(Tensor<float>& input, const Tensor<float>& bias){
    int64_t row = input.dim(0);
    int64_t col = input.dim(1);
    int64_t ld = col;
    float* input_data = input.data();
    const float* bias_data = bias.data();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int64_t r = 0; r < row; ++r){
        for(int64_t c = 0; c < col; ++c){
            input_data[r * ld + c] += bias_data[c];
        }
    }
}

void bias_naive(Tensor<double>& input, const Tensor<double>& bias){
    int64_t row = input.dim(0);
    int64_t col = input.dim(1);
    int64_t ld = col;
    double* input_data = input.data();
    const double* bias_data = bias.data();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int64_t r = 0; r < row; ++r){
        for(int64_t c = 0; c < col; ++c){
            input_data[r * ld + c] += bias_data[c];
        }
    }
}

template<typename DataType>
void Linear<DataType>::forward(const Tensor<DataType>& input, Tensor<DataType>& output){
    char transA = 'N';
    char transB = 'N';
    DataType alpha = 1.;
    DataType beta = 0.;
    int m = out_features_;
    int n = input.dim(0);
    int k = in_features_;
    DataType* A = weights_.data();
    int lda = out_features_;
    DataType* B = const_cast<DataType*>(input.data());
    int ldb = input.dim(1);
    DataType*  C = output.data();
    int ldc = output.dim(1);

    gemm(&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);

    bias_naive(output, bias_);
}

template<typename DataType>
void LinearGELU<DataType>::forward(const Tensor<DataType>& input, Tensor<DataType>& output){
    char transA = 'N';
    char transB = 'N';
    DataType alpha = 1.;
    DataType beta = 0.;
    int m = out_features_;
    int n = input.dim(0);
    int k = in_features_;
    DataType* A = weights_.data();
    int lda = out_features_;
    DataType* B = const_cast<DataType*>(input.data());
    int ldb = input.dim(1);
    DataType*  C = output.data();
    int ldc = output.dim(1);

    gemm(&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);

    bias_naive(output, bias_);
    // GELU
    gelu_exp(output.element_num(), output.data());
}

template class Tensor<float>;
template class Linear<float>;
template class LinearGELU<float>;

template class Tensor<double>;
template class Linear<double>;
template class LinearGELU<double>;