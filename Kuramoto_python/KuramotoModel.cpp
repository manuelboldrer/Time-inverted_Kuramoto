
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

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
        double x = 5 * cos(3 * t);
        double y = 5 * sin(2 * t);
        double z = 0 * sin(5 * t);
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

    return 0;
}
