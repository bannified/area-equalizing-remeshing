#pragma once

#include <string>
#include <iostream>
#include <chrono>

class ScopedTimer {
public:
    ScopedTimer(const std::string& title) : title(std::move(title)), start(std::chrono::steady_clock::now()) { }

    ~ScopedTimer() {
#if NDEBUG
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0;
        std::cout << "Scoped Timer '" << title << "' took " << elapsed << " seconds." << std::endl;
#endif
    };

private:
    std::chrono::steady_clock::time_point start;
    std::string title;
};

