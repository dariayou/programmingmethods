#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>

auto fibonacci_sequence(int a, int b) { 
    return [=]() mutable {
        int next = a + b;
        a = b;
        b = next;
        return next;
    };
}


auto random_numbers() {
    std::vector<int> numbers(11);
    std::iota(numbers.begin(), numbers.end(), 0);
    std::shuffle(numbers.begin(), numbers.end(), std::mt19937(std::random_device()()));
    
    return [numbers](int& out) mutable {
        if(numbers.empty()) return false;
        out = numbers.back();
        numbers.pop_back();
        return true;
    };
}


int main() {
    std::cout << "Задание 1: Генератор последовательностей чисел Фиббоначи" << std::endl;
    auto fibA = fibonacci_sequence(1, 1);
    std::cout << "fibA (старт 1,1):" << std::endl;
    std::cout << fibA() << std::endl;
    std::cout << fibA() << std::endl;

    auto fibB = fibonacci_sequence(0, 2);
    std::cout << "fibB (старт 0,2):" << std::endl;
    std::cout << fibB() << std::endl; 
    std::cout << fibB() << std::endl; 
    std::cout << fibB() << std::endl; 

    std::cout << "\nЗадание 2: Генератор случайных чисел" << std::endl;
    auto rnd = random_numbers();
    int num;
    int count = 1;
    
    while(true) {
        if(!rnd(num)) {
            std::cout << count << ": false" << std::endl;
            break;
        }
        std::cout << count << ": " << num << std::endl;
        count++;
    }
}
