#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <stdexcept>
#include <string>

/** Exception with customizable message. */
class GeneralException : public std::exception {
public:
    explicit GeneralException(const std::string& msg) 
        : exception(), message_(msg) {}

    virtual ~GeneralException() throw() {}

    virtual const char* what() const throw()
    {
        return message_.c_str();
    }

protected:
    const std::string message_;
};

class TreeCorrelatorException : public GeneralException {
    public:
        TreeCorrelatorException(const std::string& msg) : 
            GeneralException(msg) {}
};

class ConversionException : public GeneralException {
    public:
        ConversionException(const std::string& msg) : 
            GeneralException(msg) {}
};

#endif
