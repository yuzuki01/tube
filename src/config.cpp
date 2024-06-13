#include "tube.h"
#include "config.h"


StringList split(const String &str) {
    StringList result;
    std::istringstream iss(str);
    std::string token;

    while (iss >> token) {
        token.erase(std::remove(token.begin(), token.end(), '\r'), token.end());
        token.erase(std::remove(token.begin(), token.end(), '\n'), token.end());
        result.push_back(token);
    }

    return result;
}


Config::Config(const std::string &file_path) {
    StringList lines;
    std::ifstream file(file_path);

    {
        String line;
        if (file.is_open()) {
            while (std::getline(file, line)) {
                line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
                line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                if (!line.empty()) {
                    lines.push_back(line);
                }
            }
            file.close();
        } else {
            std::cerr << "Unable to open file: " << file_path << std::endl;
        }
    }
    /// parse lines
    for (auto &line: lines) {
        auto data = split(line);
        if (data.size() <= 1) continue;
        if (data[0] == "#") continue;
        values[data[0]] = data[1];
    }
}

template<>
int Config::get(const String &key, int _default) {
    if (values.find(key) == values.end()) return _default;
    return std::stoi(values[key]);
}

template<>
Scalar Config::get(const String &key, Scalar _default) {
    if (values.find(key) == values.end()) return _default;
    return std::stod(values[key]);
}
