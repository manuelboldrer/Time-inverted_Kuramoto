#include <iostream>
#include <vector>
#include <cmath>

class KuramotoModel
{
private:
    std::vector<std::vector<int>> circular_undirected_adjacency_matrix(int N)
    {
        // Initialize the adjacency matrix with zeros
        std::vector<std::vector<int>> adjacency_matrix(N, std::vector<int>(N, 0));

        // Fill the adjacency matrix for circular undirected edges
        for (int i = 0; i < N; ++i)
        {
            adjacency_matrix[i][(i + 1) % N] = 1;
            adjacency_matrix[(i + 1) % N][i] = 1;
        }

        return adjacency_matrix;
    }

public:
    std::vector<double> kuramoto(std::vector<double> x, double K, int N, double Omega)
    {
        // Compute the adjacency matrix
        std::vector<std::vector<int>> Ad = circular_undirected_adjacency_matrix(N);

        // Reshape x to a column vector
        std::vector<std::vector<double>> x_matrix(x.size(), std::vector<double>(1));
        for (size_t i = 0; i < x.size(); ++i)
        {
            x_matrix[i][0] = x[i];
        }

        // Compute f for the Kuramoto model
        std::vector<double> f(N);
        for (int i = 0; i < N; ++i)
        {
            double sum = 0.0;
            for (int j = 0; j < N; ++j)
            {
                sum += sin(x_matrix[j][0] * Ad[i][j] - Ad[j][i] * x_matrix[i][0]);
            }
            f[i] = Omega + (K / N) * sum;
        }

        return f;
    }
};

int main()
{
    KuramotoModel model;

    std::vector<double> x = {1.0, 2.0, 3.0, 4, 2}; // Example value for x
    double K = 0.5;                                // Example value for K
    int N = 5;                                     // Example value for N
    double Omega = 0.1;                            // Example value for Omega

    // Call the function
    std::vector<double> f = model.kuramoto(x, K, N, Omega);

    // Print the result
    std::cout << "f: ";
    for (double val : f)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}
