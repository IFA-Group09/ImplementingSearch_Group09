#include "benchmark.hpp"

#include <chrono>

Benchmark::Benchmark(std::string method, std::filesystem::path reference_path, std::filesystem::path query_path) {
   this->method = method;
   this->reference_path = reference_path;
   this->query_path = query_path;

   // setup benchmarking output file
    std::ifstream benchmark_in;
    benchmark_out.open("cpp_benchmark.csv", std::ios_base::app);

    benchmark_in.open("cpp_benchmark.csv");
    if (benchmark_in.peek() == std::ifstream::traits_type::eof()) {
	benchmark_out << "method,reference_file,reads_file,time,read_n\n";
    }
    //this->start_time = std::chrono::system_clock::now();
}

void Benchmark::write(int read_num) {
	benchmark_out << method << "," << reference_path << "," << query_path << "," << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - this->start_time).count() << "," << read_num << std::endl;
}
