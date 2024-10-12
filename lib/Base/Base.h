#pragma once
#include <cassert>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <queue>
#include <random>
#include <unordered_map>
#include <vector>

using std::deque;
using std::string;

static const bool DEBUG = false;

FILE *log_to = stderr;

namespace std {
// from boost (functional/hash):
// see http://www.boost.org/doc/libs/1_35_0/doc/html/hash/combine.html template
template <class T> inline void combine(std::size_t &seed, T const &v) {
  seed ^= hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <> struct hash<std::pair<int, int>> {
  auto operator()(const std::pair<int, int> &x) const -> size_t {
    std::size_t seed = 17;
    combine(seed, x.first);
    combine(seed, x.second);
    return seed;
  }
};

template <> struct hash<std::vector<int>> {
  auto operator()(const std::vector<int> &x) const -> size_t {
    std::size_t seed = 17;
    for (auto it : x)
      combine(seed, it);
    return seed;
  }
};
} // namespace std

template <typename T> void EraseIndex(std::vector<T> &vec, int &idx) {
  vec[idx] = vec.back();
  vec.pop_back();
  --idx;
}
