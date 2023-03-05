#include "TSP.hpp"
#include <algorithm>
#include <stack>
#include <optional>

std::ostream &operator<<(std::ostream &os, const CostMatrix &cm)
{
    for (std::size_t r = 0; r < cm.size(); ++r)
    {
        for (std::size_t c = 0; c < cm.size(); ++c)
        {
            const auto &elem = cm[r][c];
            os << (is_inf(elem) ? "INF" : std::to_string(elem)) << " ";
        }
        os << "\n";
    }
    os << std::endl;

    return os;
}

bool operator==(const tsp_solution_t &s1, const tsp_solution_t &s2)
{
    return s1.lower_bound == s2.lower_bound && s1.path == s2.path;
}

bool operator!=(const tsp_solution_t &s1, const tsp_solution_t &s2)
{
    return !(s1 == s2);
}

/* PART 1 */

/**
 * Create path from unsorted path and last 2x2 cost matrix.
 * @return The vector of consecutive vertex.
 */
path_t StageState::get_path()
{
    reduce_cost_matrix();
    NewVertex lastVertex = choose_new_vertex();
    update_cost_matrix(lastVertex.coordinates);


    NewVertex lastVertex2 = choose_new_vertex();
    update_cost_matrix(lastVertex2.coordinates);
    std::vector<size_t> final_path;

    unsorted_path_t path_sort;
    for (auto x: unsorted_path_)
    {
        path_sort.push_back(x);
    }
    path_sort.push_back(lastVertex.coordinates);
    path_sort.push_back(lastVertex2.coordinates);
    unsorted_path_.push_back(lastVertex.coordinates);
    unsorted_path_.push_back(lastVertex2.coordinates);
    std::size_t value = path_sort[0].row;

    while (final_path.size() != matrix_.size())
    {
        for (auto &x: path_sort)
        {
            if (x.row == value)
            {
                value = x.col;
                final_path.push_back(value);
            }
        }
    }

    return final_path;
}

/**
 * Get minimum values from each row and returns them.
 * @return Vector of minimum values in row.
 */
std::vector<cost_t> CostMatrix::get_min_values_in_rows() const
{
    std::vector<cost_t> minValues{};
    minValues.reserve(matrix_.size());

    for (auto &row: matrix_)
    {
        minValues.push_back(*std::min_element(row.begin(), row.end()));
    }
    return minValues;
}

/**
 * Reduce rows so that in each row at least one zero value is present.
 * @return Sum of values reduced in rows.
 */
cost_t CostMatrix::reduce_rows()
{
    std::vector<cost_t> min = get_min_values_in_rows();
    for (std::size_t i = 0; i < min.size(); i++)
    {
        for (std::size_t j = 0; j < min.size(); ++j)
        {
            if (matrix_[i][j] != INF)
                matrix_[i][j] = matrix_[i][j] - min[i];
        }
    }
    for (int &x: min)
    {
        if (x == INF)
            x = 0;
    }
    cost_t sum = std::accumulate(min.begin(), min.end(), 0);
    return sum;
}

/**
 * Get minimum values from each column and returns them.
 * @return Vector of minimum values in columns.
 */
std::vector<cost_t> CostMatrix::get_min_values_in_cols() const
{
    std::vector<cost_t> min_values;
    min_values.reserve(matrix_[0].size());

    for (std::size_t j = 0; j < std::size(matrix_[0]); ++j)
    {
        cost_t minValue = matrix_[0][j];
        for (const auto &row: matrix_)
        {
            minValue = (row[j] < minValue ? row[j] : minValue);
        }
        min_values.push_back(minValue);
    }
    return min_values;
}

/**
 * Reduces rows so that in each column at least one zero value is present.
 * @return Sum of values reduced in columns.
 */
cost_t CostMatrix::reduce_cols()
{
    std::vector<cost_t> min = get_min_values_in_cols();
    for (std::size_t i = 0; i < min.size(); i++)
    {
        for (std::size_t j = 0; j < min.size(); ++j)
        {
            if (matrix_[j][i] != INF)
                matrix_[j][i] = matrix_[j][i] - min[i];
        }
    }
    for (int &x: min)
    {
        if (x == INF)
            x = 0;
    }
    cost_t sum = std::accumulate(min.begin(), min.end(), 0);
    return sum;
}

/**
 * Get the cost of not visiting the vertex_t (@see: get_new_vertex())
 * @param row
 * @param col
 * @return The sum of minimal values in row and col, excluding the intersection value.
 */
cost_t CostMatrix::get_vertex_cost(std::size_t row, std::size_t col) const
{
    cost_t minRow = INF;
    for (std::size_t j = 0; j < matrix_[row].size(); ++j)
    {
        minRow = (matrix_[row][j] < minRow && j != col ? matrix_[row][j] : minRow);
    }
    minRow = minRow == INF ? 0 : minRow;

    cost_t minCol = INF;
    for (std::size_t i = 0; i < matrix_.size(); ++i)
    {
        minCol = (matrix_[i][col] < minCol && i != row ? matrix_[i][col] : minCol);
    }
    minCol = minCol == INF ? 0 : minCol;

    return minRow + minCol;
}

/* PART 2 */

/**
 * Choose next vertex to visit:
 * - Look for vertex_t (pair row and column) with value 0 in the current cost matrix.
 * - Get the vertex_t cost (calls get_vertex_cost()).
 * - Choose the vertex_t with maximum cost and returns it.
 * @param cm
 * @return The coordinates of the next vertex.
 */
NewVertex StageState::choose_new_vertex()
{
    NewVertex newVertex{0, -INF};
    for (std::size_t i = 0; i < matrix_.size(); ++i)
    {
        for (std::size_t j = 0; j < matrix_[i].size(); ++j)
        {
            if (matrix_[i][j] == 0)
            {
                NewVertex currVertex(vertex_t(i, j), matrix_.get_vertex_cost(i, j));
                newVertex = (currVertex.cost > newVertex.cost && currVertex.cost != INF ? currVertex : newVertex);
            }
        }
    }
    return newVertex;
}

/**
 * Update the cost matrix with the new vertex.
 * @param new_vertex
 */
void StageState::update_cost_matrix(vertex_t new_vertex)
{

    for (int &j: matrix_[new_vertex.row])
    {
        j = INF;
    }

    for (std::size_t i = 0; i < matrix_.size(); ++i)
    {
        matrix_[i][new_vertex.col] = INF;
    }


    std::size_t lastRow = INF;
    std::size_t lastCol = INF;

    for (auto v: unsorted_path_)
    {
        auto itCol = std::find_if(unsorted_path_.begin(), unsorted_path_.end(),
                                  [v](vertex_t other) { return v.col == other.row; });

        auto itRow = std::find_if(unsorted_path_.begin(), unsorted_path_.end(),
                                  [v](vertex_t other) { return v.row == other.col; });

        if (itCol == unsorted_path_.end() && itRow == unsorted_path_.end() && unsorted_path_.size() > 1)
        { continue; }

        lastCol = (itCol == unsorted_path_.end() ? v.col : lastCol);
        lastRow = (itRow == unsorted_path_.end() ? v.row : lastRow);
    }

    if (lastCol != INF && lastRow != INF && unsorted_path_.size() < matrix_.size() - 1)
    {
        matrix_[lastCol][lastRow] = INF;
    }
}

/**
 * Reduce the cost matrix.
 * @return The sum of reduced values.
 */
cost_t StageState::reduce_cost_matrix()
{
    bool need_reduce_rows = false;
    bool need_reduce_cols = false;
    for (std::size_t i = 0; i < matrix_.size(); i++)
    {
        bool row_has_zero = false;
        bool col_has_zero = false;
        for (std::size_t j = 0; j < matrix_.size(); j++)
        {
            if (matrix_[i][j] == 0)
            {
                row_has_zero = true;
            }
            if (matrix_[j][i] == 0)
            {
                col_has_zero = true;
            }
        }
        if (!row_has_zero)
        {
            need_reduce_rows = true;
        }
        if (!col_has_zero)
        {
            need_reduce_cols = true;
        }
    }

    cost_t rows_reduced = 0;
    if (need_reduce_rows)
    {
        rows_reduced = matrix_.reduce_rows();
    }
    cost_t cols_reduced = 0;
    if (need_reduce_cols)
    {
        cols_reduced = matrix_.reduce_cols();
    }

    return cols_reduced + rows_reduced;
}

/**
 * Given the optimal path, return the optimal cost.
 * @param optimal_path
 * @param m
 * @return Cost of the path.
 */
cost_t get_optimal_cost(const path_t &optimal_path, const cost_matrix_t &m)
{
    cost_t cost = 0;

    for (std::size_t idx = 1; idx < optimal_path.size(); ++idx)
    {
        cost += m[optimal_path[idx - 1]][optimal_path[idx]];
    }

    // Add the cost of returning from the last city to the initial one.
    cost += m[optimal_path[optimal_path.size() - 1]][optimal_path[0]];

    return cost;
}

/**
 * Create the right branch matrix with the chosen vertex forbidden and the new lower bound.
 * @param m
 * @param v
 * @param lb
 * @return New branch.
 */
StageState create_right_branch_matrix(cost_matrix_t m, vertex_t v, cost_t lb)
{
    CostMatrix cm(m);
    cm[v.row][v.col] = INF;
    return StageState(cm, {}, lb);
}

/**
 * Retain only optimal ones (from all possible ones).
 * @param solutions
 * @return Vector of optimal solutions.
 */
tsp_solutions_t filter_solutions(tsp_solutions_t solutions)
{
    cost_t optimal_cost = INF;
    for (const auto &s: solutions)
    {
        optimal_cost = (s.lower_bound < optimal_cost) ? s.lower_bound : optimal_cost;
    }

    tsp_solutions_t optimal_solutions;
    std::copy_if(solutions.begin(), solutions.end(),
                 std::back_inserter(optimal_solutions),
                 [&optimal_cost](const tsp_solution_t &s) { return s.lower_bound == optimal_cost; }
    );

    return optimal_solutions;
}

/**
 * Solve the TSP.
 * @param cm The cost matrix.
 * @return A list of optimal solutions.
 */
tsp_solutions_t solve_tsp(const cost_matrix_t &cm)
{

    StageState left_branch(cm);

    // The branch & bound tree.
    std::stack<StageState> tree_lifo;

    // The number of levels determines the number of steps before obtaining
    // a 2x2 matrix.
    std::size_t n_levels = cm.size() - 2;

    tree_lifo.push(left_branch);   // Use the first cost matrix as the root.

    cost_t best_lb = INF;
    tsp_solutions_t solutions;

    while (!tree_lifo.empty())
    {

        left_branch = tree_lifo.top();
        tree_lifo.pop();

        while (left_branch.get_level() != n_levels && left_branch.get_lower_bound() <= best_lb)
        {
            // Repeat until a 2x2 matrix is obtained or the lower bound is too high...

            if (left_branch.get_level() == 0)
            {
                left_branch.reset_lower_bound();
            }

            // 1. Reduce the matrix in rows and columns.
            cost_t new_cost = left_branch.reduce_cost_matrix();

            // 2. Update the lower bound and check the break condition.
            left_branch.update_lower_bound(new_cost);
            if (left_branch.get_lower_bound() > best_lb)
            {
                break;
            }

            // 3. Get new vertex and the cost of not choosing it.
            NewVertex new_vertex = left_branch.choose_new_vertex();

            // 4. Update the path - use append_to_path method.
            left_branch.append_to_path(new_vertex.coordinates);

            // 5. Update the cost matrix of the left branch.
            left_branch.update_cost_matrix(new_vertex.coordinates);

            // 6. Update the right branch and push it to the LIFO.
            cost_t new_lower_bound = left_branch.get_lower_bound() + new_vertex.cost;
            tree_lifo.push(create_right_branch_matrix(cm, new_vertex.coordinates,
                                                      new_lower_bound));
        }

        if (left_branch.get_lower_bound() <= best_lb)
        {
            // If the new solution is at least as good as the previous one,
            // save its lower bound and its path.
            best_lb = left_branch.get_lower_bound();
            path_t new_path = left_branch.get_path();
            solutions.push_back({get_optimal_cost(new_path, cm), new_path});
        }
    }

    return filter_solutions(solutions); // Filter solutions to find only optimal ones.
}