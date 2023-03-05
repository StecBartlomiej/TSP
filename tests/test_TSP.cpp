#include <gtest/gtest.h>

#include "TSP.hpp"

TEST(Test_TSP, Case_1)
{
    cost_matrix_t cm = {{INF, 10,  8,   19,  12},
                        {10,  INF, 20,  6,   3},
                        {8,   20,  INF, 4,   2},
                        {19,  6,   4,   INF, 7},
                        {12,  3,   2,   7,   INF}};
    tsp_solutions_t return_vec = solve_tsp(cm);

    tsp_solutions_t solution = {{32, {2, 3, 4, 1, 0}},
                                {32, {1, 4, 3, 2, 0}}};
    EXPECT_EQ(return_vec, solution);
}

TEST(Test_TSP, Case_2)
{
    cost_matrix_t cm{
            {INF, 12,  3,   45,  6},
            {78,  INF, 90,  21,  3},
            {5,   56,  INF, 23,  98},
            {12,  6,   8,   INF, 34},
            {3,   98,  3,   2,   INF}
    };
    tsp_solutions_t return_vec = solve_tsp(cm);

    tsp_solutions_t solution = {{30, {4, 3, 2, 0, 1}}};
    EXPECT_EQ(return_vec, solution);
}

TEST(Test_TSP, Case_3)
{
    cost_matrix_t cm{
            {INF,  3,  4,  2,  7},
            {3,  INF,  4,  6,  3},
            {4,  4,  INF,  5,  8},
            {2,  6,  5,  INF,  6},
            {7,  3,  8,  6,  INF},
    };
    tsp_solutions_t return_vec = solve_tsp(cm);

    tsp_solutions_t solution = {{19, {4, 3, 0, 2, 1}}, {19, {1, 2, 0, 3, 4}}};
    EXPECT_EQ(return_vec, solution);
}