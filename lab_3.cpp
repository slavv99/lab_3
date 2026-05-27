#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <cstddef>
#include <string>

template <typename T>
class Matrix : public std::vector<std::vector<T>> {
protected:
    std::size_t total_r;
    std::size_t total_c;
public:
    explicit Matrix(std::size_t r, std::size_t c)
    : std::vector<std::vector<T>>(r, std::vector<T>(c, static_cast<T>(0))), total_r(r), total_c(c) {
        if (r==0 || c==0){
            throw std::invalid_argument("Matrix dimensions must be positive.");
        }
    }

    Matrix(const Matrix&) = default;
    Matrix(Matrix&&) noexcept = default;
    Matrix& operator = (const Matrix&) = default;
    Matrix& operator = (Matrix&&) noexcept = default;
    virtual ~Matrix() = default;

    std::size_t rows_qty() const noexcept {return total_r;}
    std::size_t cols_qty() const noexcept {return total_c;}


    // элементарные операции метода Гаусса
    void row_swap(std::size_t idx1, std::size_t idx2) noexcept {
        std::swap(this->at(idx1), this->at(idx2));
    }

    // вычитание строк для прямого хода (зануление снизу под опорным элементом)
    void eliminate_under(std::size_t target_r, std::size_t pivot_r, std::size_t col_idx) noexcept {
        if (target_r > pivot_r) {
            T factor = this->at(target_r).at(col_idx) / this->at(pivot_r).at(col_idx);
            for (std::size_t c = col_idx; c < total_c; ++c) {
                this->at(target_r).at(c) -= this->at(pivot_r).at(c) * factor;
            }
        }
    }

    // выбор ведущего элемента по столбцу для исключения погрешностей
    std::size_t find_pivot_row(std::size_t target_c, std::size_t start_r) const noexcept {
        std::size_t lead_r = start_r;
        for (std::size_t i = start_r; i < total_r; ++i) {
            if (std::abs(this->at(lead_r).at(target_c)) < std::abs(this->at(i).at(target_c))) {
                lead_r = i;
            }
        }
        return lead_r;
    }
};

// производный класс слау, реализующий двухфазный Метод Гаусса
template <typename T>
class SystemSolver : public Matrix<T> {
private:
    std::size_t vars_num; // Число переменных (столбцов в матрице)
    const T EPS = static_cast<T>(1e-9);

    // прямой ход (приведение к ступенчатому виду)
    std::size_t forward_elimination() noexcept {
        std::size_t pivot_row = 0;
        for(std::size_t c = 0; c < vars_num && pivot_row < this->total_r; ++c) {
            std::size_t lead_r = this->find_pivot_row(c, pivot_row);
            if (std::abs(this->at(lead_r).at(c)) < EPS) {
                continue;
            }
            this->row_swap(lead_r, pivot_row);
            for (std::size_t r = pivot_row + 1; r < this->total_r; ++r) {
                this->eliminate_under(r, pivot_row, c);
            }
            ++pivot_row;
        }
        return pivot_row;
    }

    // обратный ход (нахождение решения)
    void back_substitution(std::size_t rank) {
        // определение статуса переменных (базисные или свободные)
        std::vector<bool> is_free(vars_num, true);
        std::vector<int> pivot_col_at_row(this->total_r, -1);

        for(std::size_t i = 0; i < rank; ++i) {
            for(std::size_t j = 0; j < vars_num; ++j) {
                if (std::abs(this->at(i).at(j)) > EPS) {
                    is_free[j] = false;
                    pivot_col_at_row[i] = static_cast<int>(j);
                    break;
                }
            }
        }
        if (rank == vars_num) {
            // 1 случай) единственное решение (проход снизу вверх)
            std::vector<T> solution(vars_num, static_cast<T>(0));

            for (int i = static_cast<int>(rank) - 1; i >= 0; --i) {
                std::size_t p_col = static_cast<size_t>(pivot_col_at_row[i]);
                T sum = this->at(i).at(vars_num);

                for (std::size_t j = p_col + 1; j < vars_num; ++j) {
                    sum -= this->at(i).at(j) * solution[j];
                }
                solution[p_col] = sum / this->at(i).at(p_col);
            }

            std::cout << "\n[Результат]: Единственное решение системы (Стандартный обратный ход):\n";
            for (std::size_t i = 0; i < vars_num; ++i) {
                std::cout << "x" << i + 1 << " = " << solution[i] << std::endl;
        
            }
        }
        else {
            // 2 случай) бесконечно много решений (параметрический обратный ход)
            std::cout << "\n[Результат]: Бесконечно много решений (Стандартный обратный ход): \n";

            //идем снизу вверх по ступеням матрицы выражения зависимостей
            for (int i = static_cast<int>(rank) - 1; i >= 0; --i) {
                std::size_t p_col = static_cast<size_t>(pivot_col_at_row[i]);
                T lead_coeff = this->at(i).at(p_col);
                std::cout << "x" << p_col + 1 << " = " << (this->at(i).at(vars_num) / lead_coeff);

                // выводим влияние всех последующих переменных (и базисных, и свободных)
                for(std::size_t j = p_col + 1; j < vars_num; ++j) {
                    T current_val = this->at(i).at(j);
                    if (std::abs(current_val) > EPS) {
                        T final_coeff = -current_val / lead_coeff;
                        char sign = (final_coeff > 0) ? '+' : '-';

                        std::cout << " " << sign << " " << std::abs(final_coeff);
                        std::cout << (is_free[j] ? "*t" : "*x") << j + 1;
                    }
                }
                std::cout << std::endl;
            }
            for (std::size_t j = 0; j < vars_num; ++j) {
                if (is_free[j]) {
                    std::cout << "x" << j + 1 << " = t" << j + 1 << " (свободный параметр) " << std::endl;

                }
            }

        }
    }
    
public:
    explicit SystemSolver(const Matrix<T>& matrix_A, const Matrix<T>& vector_b)
        : Matrix<T>(matrix_A.rows_qty(), matrix_A.cols_qty() + 1), vars_num(matrix_A.cols_qty()) {
        
        if (matrix_A.rows_qty() != vector_b.rows_qty() || vector_b.cols_qty() != 1) {
            throw std::invalid_argument("Mismatched SLAU components.");
        }

        for (std::size_t i = 0; i < matrix_A.rows_qty(); ++i) {
            for(std::size_t j = 0; j < matrix_A.cols_qty(); ++j) {
                this->at(i).at(j) = matrix_A.at(i).at(j);
            }
            this->at(i).at(vars_num) = vector_b.at(i).at(0);
        }
    }
    void print_system() const noexcept {
        for (const auto& row : *this) {
            std::size_t c_pos = 0;
            for (const auto& val : row) {
                if (c_pos == vars_num) std::cout << "| ";
                std::cout << val << "\t";
                c_pos++;
            }
            std::cout << "\n";
        }
    }
    // управляющая функция
    void execute_gauss() {
        std::size_t rank = forward_elimination();
        // проверка совместимости (Теорема Кронекера-Капелли)
        for (std::size_t i = rank; i < this->total_r; ++i) {
            if (std::abs(this->at(i).at(vars_num)) > EPS) {
                throw std::runtime_error("Incompatible system of equations (No solutions).");
            }
        }
        //
        back_substitution(rank);
    }
};
int main() {
    setlocale(LC_ALL, "Russian");
    std::size_t r_count = 0, c_count = 0;
    std::cout << "Введите число строк: ";
    if (!(std::cin >> r_count)) return 1;
    std::cout << "Введите число столбцов: ";
    if (!(std::cin >> c_count)) return 1;

    try {
        Matrix<double> Z(r_count, c_count);
        Matrix<double> b(r_count, 1);

        std::cout << "Введите коэффициенты матрицы :\n";
        for (std::size_t i = 0; i < r_count; ++i) {
            for(std::size_t j = 0; j < c_count; ++j) {
                std::cout << "Z[" << i << "][" << j << "] = ";
                if (!(std::cin >> Z.at(i).at(j))) 
                    throw std::invalid_argument("Invalid data format.");
            }
        }

        std::cout << "Введите свободные члены :\n";
        for (std::size_t i = 0; i < r_count; ++i) {
            std::cout << "b[" << i << "][0] = ";
            if (!(std::cin >> b.at(i).at(0)))
                throw std::invalid_argument("Invalid data format.");
        }

        SystemSolver<double> solver(Z, b);

        std::cout << "\nРасширенная матрица системы: \n";
        solver.print_system();
        solver.execute_gauss();
    }
    catch (const std::bad_alloc&) {
        std::cerr <<"[Memory Alloc Error] : Out of memory space.\n";
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "\n[Argument Error]: " << e.what() << "\n";
    }
    catch (const std::runtime_error& e) {
        std::cerr << "\n[Runtime Math Error]: " <<e.what() << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "\n[Standard Exception]: " << e.what() << "\n";
    }
    return 0;
}
