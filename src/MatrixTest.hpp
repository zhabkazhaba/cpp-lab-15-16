#ifndef CPP_LAB_15_16_MATRIXTEST_HPP
#define CPP_LAB_15_16_MATRIXTEST_HPP

#include "Matrix.hpp"

enum TEST_MODE {
    SUM,
    MUL
};

class MatrixTest {
private:
    using testResult = std::tuple<std::vector<int>, std::vector<long long>, std::vector<long long>>; // Variable, serial vector, parallel vector
    using testRange = std::pair<unsigned int, unsigned int>;
public:
    static testResult runConstThreads(TEST_MODE mode, testRange range, unsigned int step,
                               unsigned int num_threads = std::thread::hardware_concurrency()) {
        std::vector<int> vars;
        std::vector<long long> serial_results, parallel_results;
        for (unsigned int i = range.first; i < range.second + 1; i += step) {
            Matrix<float> Mat1(i, i);
            Mat1.initializeWithNumbers();
            Matrix<float> Mat2(i, i);
            Mat2.initializeWithNumbers();
            Matrix<float> Mat3(i, i);
            Mat3.initializeWithNumbers();
            auto start = std::chrono::high_resolution_clock::now();
            if (mode == TEST_MODE::SUM)
                Mat1.parallelSum(Mat2, num_threads);
            else
                Mat1.parallelMul(Mat2, num_threads);
            auto end = std::chrono::high_resolution_clock::now();
            long long t1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            start = std::chrono::high_resolution_clock::now();
            if (mode == TEST_MODE::SUM)
                Mat2.serialSum(Mat3);
            else
                Mat2.serialMul(Mat3);
            end = std::chrono::high_resolution_clock::now();
            long long t2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            vars.emplace_back(i);
            serial_results.emplace_back(t2);
            parallel_results.emplace_back(t1);
            std::cout << "Checked i = " << i << std::endl;
        }
        testResult res = std::make_tuple(vars, serial_results, parallel_results);
        return res;
    }

    static testResult runConstSize(TEST_MODE mode, testRange range, unsigned int step, unsigned int matrix_size = 100) {
        std::vector<int> vars;
        std::vector<long long> serial_results, parallel_results;
        for (unsigned int i = range.first; i < range.second + 1; i += step) {
            Matrix<float> Mat1(matrix_size, matrix_size);
            Mat1.initializeWithNumbers();
            Matrix<float> Mat2(matrix_size, matrix_size);
            Mat2.initializeWithNumbers();
            Matrix<float> Mat3(matrix_size, matrix_size);
            Mat3.initializeWithNumbers();
            auto start = std::chrono::high_resolution_clock::now();
            if (mode == TEST_MODE::SUM)
                Mat1.parallelSum(Mat2, i);
            else
                Mat1.parallelMul(Mat2, i);
            auto end = std::chrono::high_resolution_clock::now();
            long long t1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            start = std::chrono::high_resolution_clock::now();
            if (mode == TEST_MODE::SUM)
                Mat2.serialSum(Mat3);
            else
                Mat2.serialMul(Mat3);
            end = std::chrono::high_resolution_clock::now();
            long long t2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            vars.emplace_back(i);
            serial_results.emplace_back(t2);
            parallel_results.emplace_back(t1);
            std::cout << "Checked i = " << i << std::endl;
        }
        testResult res = std::make_tuple(vars, serial_results, parallel_results);
        return res;
    }

    static void printResults(const testResult& result) {
        if (std::get<0>(result).empty()) {
            std::cerr << "Result vector is empty." << std::endl;
            return;
        }
        std::cout << "Var:  " << "Serial: " << "Parallel: " << "Speedup: " << std::endl;
        for (unsigned int i = 0; i < std::get<0>(result).size(); i++) {
            int var = std::get<0>(result)[i];
            unsigned int t1 = std::get<1>(result)[i];
            unsigned int t2 = std::get<2>(result)[i];
            std::cout << var << " " << t1 << " " << t2 << " " << float(t1) / float(t2) << std::endl;
        }
    }
};

#endif //CPP_LAB_15_16_MATRIXTEST_HPP
