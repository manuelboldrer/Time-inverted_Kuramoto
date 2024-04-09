#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <limits>

class KuramotoModel
{

public:
    double kuramoto(std::vector<double> theta, double K, double Omega)
    {
        // Compute f for the Kuramoto model
        double f;
        f = Omega - K * (sin(theta[1] - theta[0]) + sin(theta[2] - theta[0]));
        return f;
    }

    std::vector<double> fun(double X, double alpha)
    {

        X = alpha * X;
        double t = fmod(X, 2 * M_PI);
        double x = 35 * cos(7 * t);
        double y = 45 * sin(6 * t);
        double z = 5 * sin(5 * t);
        std::vector<double> result = {x, y, z};
        return result;
    }

    // Function to calculate the distance between two points in 3D space
    double distance(double x1, double y1, double z1, double x2, double y2, double z2)
    {
        return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
    }

    // Function to find the closest point on the curve to the given point
    // Function to find the closest point on the curve to the given point
    void closestPointOnCurve(double x0, double y0, double z0, double alpha, double &closestX, double &closestY, double &closestZ, double &closestT)
    {
        const int numPoints = 10000;                         // Number of points to sample along the curve
        double minDist = std::numeric_limits<double>::max(); // Initialize minimum distance to a large value
        double prevClosestT = closestT;                      // Store the previous closest parameter value
        closestT = 0.0;

        // Iterate through parameter values to sample points on the curve
        for (int i = 0; i <= numPoints; ++i)
        {
            double t = i / static_cast<double>(numPoints);  // Adjust t from 0 to 1
            std::vector<double> curvePoint = fun(t, alpha); // alpha value
            double x_t = curvePoint[0];
            double y_t = curvePoint[1];
            double z_t = curvePoint[2];
            double dist = distance(x0, y0, z0, x_t, y_t, z_t);

            // Update closest point if the distance is smaller
            if (dist < minDist)
            {
                minDist = dist;
                closestT = t;
                closestX = x_t;
                closestY = y_t;
                closestZ = z_t;
            }
        }

        // Check if the new closest point has a closestT within a small distance of the previous closestT
        if (std::abs(closestT - prevClosestT) < 1.0)
        {
            // If so, use the new closest point
            // No action needed, as the current closest point is already selected
        }
        else
        {
            // If not, revert to the previous closest point
            closestT = prevClosestT;
            std::vector<double> curvePoint = fun(prevClosestT, alpha);
            closestX = curvePoint[0];
            closestY = curvePoint[1];
            closestZ = curvePoint[2];
        }
    }
};

int main()
{
    KuramotoModel kuramoto;
    // TODO: convert position to theta.
    double givenX = -30.85311913666038; // Example: Given x-coordinate of the point
    double givenY = 47.91416182741961;  // Example: Given y-coordinate of the point
    double givenZ = 3.586899449499661;  // Example: Given z-coordinate of the point
    double alpha = 4.0;                 // Example: Alpha parameter for the

    double closestX, closestY, closestZ, closestT;
    kuramoto.closestPointOnCurve(givenX, givenY, givenZ, alpha, closestX, closestY, closestZ, closestT);

    std::cout << "Closest point on the curve: (" << closestX << ", " << closestY << ", " << closestZ << ")" << std::endl;
    std::cout << "Parameter value for the closest point: " << closestT << std::endl;

    std::vector<double> theta = {closestT, 2.09, 4.18}; // Example value for theta
    double K = 0.5;                                     // Example value for K
    double Omega = 0;                                   // Example value for Omega
    double dt = 0.1;                                    // Example value for dt

    // Call the function kuramoto
    double f = kuramoto.kuramoto(theta, K, Omega);

    // Print the result
    std::cout << "f: " << f << std::endl;

    // Calculate the next theta value
    double thetanext_i = theta[0] + dt * f;

    // Calculate the next x, y, z values
    // double alpha = 1.0; // Example value for alpha
    std::vector<double> next_values = kuramoto.fun(theta[0], alpha);
    double xnext_i = next_values[0];
    double ynext_i = next_values[1];
    double znext_i = next_values[2];

    // Print the results
    std::cout << "thetanext_i: " << thetanext_i << std::endl;
    std::cout << "xnext_i: " << xnext_i << std::endl;
    std::cout << "ynext_i: " << ynext_i << std::endl;
    std::cout << "znext_i: " << znext_i << std::endl;

    return 0;
}