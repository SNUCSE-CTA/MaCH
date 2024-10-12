#pragma once
#include "Base/Base.h"

deque<string> parse(string line, const string &del) {
    deque<string> ret;

    size_t pos = 0;
    string token;
    while ((pos = line.find(del)) != string::npos) {
        token = line.substr(0, pos);
        ret.push_back(token);
        line.erase(0, pos + del.length());
    }
    ret.push_back(line);
    return ret;
}

static std::streampos fileSize(const char *filePath) {
    std::streampos fsize = 0;
    std::ifstream file(filePath, std::ios::binary);

    fsize = file.tellg();
    file.seekg(0, std::ios::end);
    fsize = file.tellg() - fsize;
    file.close();

    return fsize;
}

bool file_ok(std::string &filepath) {
    std::ifstream infile(filepath);
    return infile.good();
}

bool CreateDirectory(const std::string &dirName) {
    std::error_code err;
    if (!std::filesystem::create_directories(dirName, err)) {
        if (std::filesystem::exists(dirName)) {
            return true;
        }
        printf("CREATEDIR: FAILED to create [%s], err:%s\n", dirName.c_str(),
               err.message().c_str());
        return false;
    }
    return true;
}

template <typename TType>
void LogVector(const std::vector<TType> &vec, std::ostream &out = std::cout,
               std::string name = "", bool endline = true) {
    typename std::vector<TType>::const_iterator it;
    if (!name.empty()) out << "[" << name << "]: ";
    out << "(";
    for (it = vec.begin(); it != vec.end(); ++it) {
        if (it != vec.begin()) out << ",";
        out << (*it);
    }
    out << ")";
    if (endline) out << "\n";
}

template <typename TType>
void PrintVector(const std::vector<TType> &vec, std::ostream &out = std::cout,
                 bool binary = true, bool print_num = true) {
    typename std::vector<TType>::const_iterator it;
    int num_elements = vec.size();
    if (print_num)
        out.write(reinterpret_cast<const char *>(&num_elements),
                  sizeof(num_elements));
    if (!vec.empty()) {
        out.write(reinterpret_cast<const char *>(vec.data()),
                  num_elements * sizeof(TType));
    }
}

template <typename TType>
bool ReadVector(std::vector<TType> &v, std::ifstream &inFile,
                bool binary = true, int num_elements = -1) {
    if (binary) {
        if (num_elements < 0) {
            if (!inFile.read(reinterpret_cast<char *>(&num_elements),
                             sizeof(num_elements))) {
                return false;
            }
        }
        v.resize(num_elements);
        if (num_elements > 0) {
            inFile.read(reinterpret_cast<char *>(v.data()),
                        num_elements * sizeof(TType));
        }
        return true;
    } else {
        fprintf(stderr, "ReadVector: not implemented\n");
        exit(EXIT_FAILURE);
    }
    return true;
}
