#ifndef TUBE_CONFIG_H
#define TUBE_CONFIG_H

typedef std::string String;
typedef std::vector<String> StringList;

StringList split(const String &str);

class Config {
protected:
    std::unordered_map<std::string, std::string> values;
public:
    explicit Config(const std::string &file_path);

    template<class DataType>
    DataType get(const String &key, DataType _default);
};

#endif //TUBE_CONFIG_H
