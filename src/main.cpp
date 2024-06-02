#include <iostream>
#include <chrono>
#include "Matrix.hpp"
#include "MatrixTest.hpp"

using mt = MatrixTest;

int main() {
    std::cout << std::boolalpha;
    std::pair<unsigned int, unsigned int> range1 = {1,700};
    std::pair<unsigned int, unsigned int> range2 = {1, 7};

    //mt::printResults(mt::runConstThreads(TEST_MODE::SUM, range1, 50, 4));
    mt::printResults(mt::runConstThreads(TEST_MODE::MUL, range1, 20, 6));
    std::cout << "--------------------------------" << std::endl;
    //mt::printResults(mt::runConstThreads(TEST_MODE::MUL, range, 5, 4));
    return 0;
}
