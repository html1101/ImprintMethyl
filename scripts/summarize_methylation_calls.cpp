#include <array>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

namespace {

struct ContextCounts {
    std::uint64_t methylated = 0;
    std::uint64_t unmethylated = 0;
    std::uint64_t covered_sites = 0;
};

std::string fraction(std::uint64_t methylated, std::uint64_t total, int precision) {
    if (total == 0) {
        return "NA";
    }
    char buffer[64];
    std::snprintf(
        buffer, sizeof(buffer), precision == 6 ? "%.6g" : "%.8g",
        static_cast<double>(methylated) / static_cast<double>(total)
    );
    return buffer;
}

bool parse_count(std::string_view value, std::uint64_t& count) {
    if (value.empty()) {
        return false;
    }
    const auto result = std::from_chars(value.data(), value.data() + value.size(), count);
    return result.ec == std::errc{} && result.ptr == value.data() + value.size();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: summarize_methylation_calls CROSS STRAND MIN_COVERAGE SUMMARY_PATH\n";
        return 2;
    }

    const std::string cross = argv[1];
    const std::string strand = argv[2];
    const std::uint64_t min_coverage = std::stoull(argv[3]);
    const std::string summary_path = argv[4];
    if (cross != "F1i" && cross != "F1r") {
        std::cerr << "unsupported cross label: " << cross << '\n';
        return 2;
    }

    const std::string expected_header =
        "chr\tpos\ttype\t" + cross + "_129_met\t" + cross + "_129_unmet\t" +
        cross + "_Cast_met\t" + cross + "_Cast_unmet\t" +
        cross + "_unassigned_met\t" + cross + "_unassigned_unmet";
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::string line;
    if (!std::getline(std::cin, line) || line != expected_header) {
        std::cerr << "unexpected methylation-call table header\n";
        return 2;
    }

    std::array<std::array<ContextCounts, 3>, 3> summary{};
    const std::array<std::string_view, 3> contexts = {"CG", "CHG", "CHH"};
    const std::array<std::string_view, 3> alleles = {"129", "Cast", "unassigned"};
    std::uint64_t rows = 0;
    std::uint64_t written_sites = 0;

    while (std::getline(std::cin, line)) {
        ++rows;
        std::array<std::string_view, 9> fields;
        std::size_t start = 0;
        bool valid_columns = true;
        for (std::size_t index = 0; index < fields.size(); ++index) {
            const auto tab = line.find('\t', start);
            if (index + 1 == fields.size()) {
                if (tab != std::string::npos) {
                    valid_columns = false;
                    break;
                }
                fields[index] = std::string_view(line).substr(start);
            } else {
                if (tab == std::string::npos) {
                    valid_columns = false;
                    break;
                }
                fields[index] = std::string_view(line).substr(start, tab - start);
                start = tab + 1;
            }
        }
        if (!valid_columns) {
            std::cerr << "record " << rows << ": expected 9 columns\n";
            return 2;
        }

        std::array<std::uint64_t, 6> values{};
        for (std::size_t index = 0; index < values.size(); ++index) {
            if (!parse_count(fields[index + 3], values[index])) {
                std::cerr << "record " << rows << ": invalid count\n";
                return 2;
            }
        }

        std::size_t context_index = contexts.size();
        for (std::size_t index = 0; index < contexts.size(); ++index) {
            if (fields[2] == contexts[index]) {
                context_index = index;
                break;
            }
        }
        if (context_index == contexts.size()) {
            continue;
        }

        for (std::size_t allele_index = 0; allele_index < alleles.size(); ++allele_index) {
            const auto methylated = values[2 * allele_index];
            const auto unmethylated = values[2 * allele_index + 1];
            auto& counts = summary[context_index][allele_index];
            counts.methylated += methylated;
            counts.unmethylated += unmethylated;
            if (methylated + unmethylated > 0) {
                ++counts.covered_sites;
            }
        }

        if (context_index != 0) {
            continue;
        }
        const auto m129 = values[0];
        const auto u129 = values[1];
        const auto mcast = values[2];
        const auto ucast = values[3];
        const auto cov129 = m129 + u129;
        const auto covcast = mcast + ucast;
        if (cov129 + covcast < min_coverage) {
            continue;
        }

        const auto maternal_m = cross == "F1i" ? m129 : mcast;
        const auto maternal_u = cross == "F1i" ? u129 : ucast;
        const auto paternal_m = cross == "F1i" ? mcast : m129;
        const auto paternal_u = cross == "F1i" ? ucast : u129;
        const auto maternal_coverage = maternal_m + maternal_u;
        const auto paternal_coverage = paternal_m + paternal_u;

        std::cout << cross << '\t' << strand << '\t' << fields[0] << '\t' << fields[1]
                  << '\t' << fields[2] << '\t' << m129 << '\t' << u129 << '\t'
                  << cov129 << '\t' << fraction(m129, cov129, 6) << '\t'
                  << mcast << '\t' << ucast << '\t' << covcast << '\t'
                  << fraction(mcast, covcast, 6) << '\t' << maternal_m << '\t'
                  << maternal_u << '\t' << maternal_coverage << '\t'
                  << fraction(maternal_m, maternal_coverage, 6) << '\t'
                  << paternal_m << '\t' << paternal_u << '\t' << paternal_coverage
                  << '\t' << fraction(paternal_m, paternal_coverage, 6) << '\t'
                  << values[4] << '\t' << values[5] << '\n';
        ++written_sites;
    }

    std::ofstream summary_output(summary_path);
    if (!summary_output) {
        std::cerr << "cannot write summary file: " << summary_path << '\n';
        return 2;
    }
    for (std::size_t context_index = 0; context_index < contexts.size(); ++context_index) {
        for (std::size_t allele_index = 0; allele_index < alleles.size(); ++allele_index) {
            const auto& counts = summary[context_index][allele_index];
            const auto total = counts.methylated + counts.unmethylated;
            summary_output << cross << '\t' << strand << '\t' << contexts[context_index]
                           << '\t' << alleles[allele_index] << '\t' << counts.methylated
                           << '\t' << counts.unmethylated << '\t' << counts.covered_sites
                           << '\t' << fraction(counts.methylated, total, 8) << '\n';
        }
    }
    if (!summary_output) {
        std::cerr << "failed while writing summary file: " << summary_path << '\n';
        return 2;
    }
    std::cerr << cross << '/' << strand << ": " << written_sites
              << " covered CpG rows written from " << rows << " records\n";
    return std::cout ? 0 : 2;
}
