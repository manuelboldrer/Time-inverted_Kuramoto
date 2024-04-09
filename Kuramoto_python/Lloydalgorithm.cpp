#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

std::vector<std::pair<double, double>> points_inside_circle(std::pair<double, double> robot_pos, double radius, double step_size)
{
    double x_center = robot_pos.first;
    double y_center = robot_pos.second;
    int x_min = static_cast<int>((x_center - radius) / step_size);
    int x_max = static_cast<int>((x_center + radius) / step_size);
    int y_min = static_cast<int>((y_center - radius) / step_size);
    int y_max = static_cast<int>((y_center + radius) / step_size);

    std::vector<double> x_coords, y_coords;
    for (int i = x_min; i <= x_max; ++i)
        x_coords.push_back(i * step_size);
    for (int j = y_min; j <= y_max; ++j)
        y_coords.push_back(j * step_size);

    std::vector<std::pair<double, double>> points;
    for (auto x : x_coords)
    {
        for (auto y : y_coords)
        {
            double distance = std::sqrt(std::pow((x - x_center), 2) + std::pow((y - y_center), 2));
            if (distance <= radius)
                points.push_back(std::make_pair(x, y));
        }
    }
    return points;
}

std::vector<std::pair<double, double>> find_closest_points(const std::pair<double, double> &robot_pos, const std::vector<std::pair<double, double>> &points, const std::vector<std::pair<double, double>> &neighbors)
{
    std::vector<double> distances_to_robot;
    for (const auto &point : points)
    {
        double distance = std::sqrt(std::pow((point.first - robot_pos.first), 2) + std::pow((point.second - robot_pos.second), 2));
        distances_to_robot.push_back(distance);
    }

    std::vector<std::vector<double>> distances_to_neighbors;
    for (const auto &point : points)
    {
        std::vector<double> distances;
        for (const auto &neighbor : neighbors)
        {
            double distance = std::sqrt(std::pow((point.first - neighbor.first), 2) + std::pow((point.second - neighbor.second), 2));
            distances.push_back(distance);
        }
        distances_to_neighbors.push_back(distances);
    }

    std::vector<std::pair<double, double>> closer_points;
    for (size_t i = 0; i < points.size(); ++i)
    {
        bool is_closer = true;
        for (size_t j = 0; j < neighbors.size(); ++j)
        {
            if (distances_to_robot[i] >= distances_to_neighbors[i][j])
            {
                is_closer = false;
                break;
            }
        }
        if (is_closer)
            closer_points.push_back(points[i]);
    }

    return closer_points;
}

std::vector<double> compute_scalar_value(const std::vector<double> &x_test, const std::vector<double> &y_test, const std::pair<double, double> &destination, double beta)
{
    std::vector<double> scalar_values;
    for (size_t i = 0; i < x_test.size(); ++i)
    {
        double distance = std::sqrt(std::pow((x_test[i] - destination.first), 2) + std::pow((y_test[i] - destination.second), 2));
        double scalar_value = std::exp(-distance / beta);
        scalar_values.push_back(scalar_value);
    }
    return scalar_values;
}

std::vector<std::pair<double, double>> account_encumbrance(const std::vector<std::pair<double, double>> &points, const std::pair<double, double> &robot_pos, const std::vector<std::pair<double, double>> &neighbors, const std::vector<double> &size_neighbors, double encumbrance)
{
    std::vector<size_t> index;
    double robot_x = robot_pos.first;
    double robot_y = robot_pos.second;

    for (size_t j = 0; j < neighbors.size(); ++j)
    {
        double delta_x = robot_x - neighbors[j].first;
        double delta_y = robot_y - neighbors[j].second;

        if (std::abs(delta_y) < 0.001)
        {
            delta_y = 0.001;
        }
        if (std::abs(delta_x) < 0.001)
        {
            delta_x = 0.001;
        }

        double m = delta_y / delta_x;
        if (std::abs(m) < 0.001)
        {
            m = 0.001;
        }

        double xm = 0.5 * (robot_x + neighbors[j].first);
        double ym = 0.5 * (robot_y + neighbors[j].second);
        double dm = std::sqrt(std::pow(xm - robot_x, 2) + std::pow(ym - robot_y, 2));

        if (dm < size_neighbors[j] + encumbrance)
        {
            std::vector<double> uvec = {delta_x, delta_y};
            double norm_uvec = std::sqrt(std::pow(delta_x, 2) + std::pow(delta_y, 2));
            uvec[0] /= norm_uvec;
            uvec[1] /= norm_uvec;

            double solx = xm + (size_neighbors[j] + encumbrance - dm) * uvec[0];
            double soly = ym + (size_neighbors[j] + encumbrance - dm) * uvec[1];

            if (robot_y + (1 / m) * (robot_x - solx) - soly > 0)
            {
                for (size_t i = 0; i < points.size(); ++i)
                {
                    if (points[i].second + (1 / m) * (points[i].first - solx) - soly < 0)
                    {
                        index.push_back(i);
                    }
                }
            }
            else
            {
                for (size_t i = 0; i < points.size(); ++i)
                {
                    if (points[i].second + (1 / m) * (points[i].first - solx) - soly > 0)
                    {
                        index.push_back(i);
                    }
                }
            }
        }
    }

    std::vector<std::pair<double, double>> new_points;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (std::find(index.begin(), index.end(), i) == index.end())
        {
            new_points.push_back(points[i]);
        }
    }
    return new_points;
}
std::pair<std::pair<double, double>, std::pair<double, double>> get_centroid(
    std::pair<double, double> robot_pos,
    double radius,
    double step_size,
    const std::vector<std::pair<double, double>> &neighbors,
    const std::vector<double> &size_neighbors,
    double encumbrance,
    const std::pair<double, double> &destination,
    double beta)
{
    // Get points inside the circle
    std::vector<std::pair<double, double>> circle_points = points_inside_circle(robot_pos, radius, step_size);

    std::vector<std::pair<double, double>> voronoi_circle_intersection;
    if (!neighbors.empty())
    {
        // Compute the Voronoi cell
        voronoi_circle_intersection = find_closest_points(robot_pos, circle_points, neighbors);
        // Account encumbrance
        voronoi_circle_intersection = account_encumbrance(voronoi_circle_intersection, robot_pos, neighbors, size_neighbors, encumbrance);
        if (voronoi_circle_intersection.empty())
        {
            voronoi_circle_intersection.push_back(robot_pos);
        }
    }
    else
    {
        voronoi_circle_intersection = circle_points;
    }

    std::vector<double> x_in, y_in;
    for (const auto &point : voronoi_circle_intersection)
    {
        x_in.push_back(point.first);
        y_in.push_back(point.second);
    }

    std::vector<double> x_in_no_neigh, y_in_no_neigh;
    for (const auto &point : circle_points)
    {
        x_in_no_neigh.push_back(point.first);
        y_in_no_neigh.push_back(point.second);
    }

    // Compute scalar values
    std::vector<double> scalar_values = compute_scalar_value(x_in, y_in, destination, beta);
    std::vector<double> scalar_values_no_neigh = compute_scalar_value(x_in_no_neigh, y_in_no_neigh, destination, beta);

    // Compute the weighted centroid
    double sum_x_in_times_scalar_values = 0.0;
    double sum_y_in_times_scalar_values = 0.0;
    double sum_scalar_values = 0.0;

    for (size_t i = 0; i < x_in.size(); ++i)
    {
        sum_x_in_times_scalar_values += x_in[i] * scalar_values[i];
        sum_y_in_times_scalar_values += y_in[i] * scalar_values[i];
        sum_scalar_values += scalar_values[i];
    }

    std::pair<double, double> centroid = std::make_pair(sum_x_in_times_scalar_values / sum_scalar_values, sum_y_in_times_scalar_values / sum_scalar_values);

    // Compute the centroid without neighbors
    double sum_x_in_no_neigh_times_scalar_values = 0.0;
    double sum_y_in_no_neigh_times_scalar_values = 0.0;
    double sum_scalar_values_no_neigh = 0.0;

    for (size_t i = 0; i < x_in_no_neigh.size(); ++i)
    {
        sum_x_in_no_neigh_times_scalar_values += x_in_no_neigh[i] * scalar_values_no_neigh[i];
        sum_y_in_no_neigh_times_scalar_values += y_in_no_neigh[i] * scalar_values_no_neigh[i];
        sum_scalar_values_no_neigh += scalar_values_no_neigh[i];
    }

    std::pair<double, double> centroid_no_neighbors = std::make_pair(sum_x_in_no_neigh_times_scalar_values / sum_scalar_values_no_neigh, sum_y_in_no_neigh_times_scalar_values / sum_scalar_values_no_neigh);

    return std::make_pair(centroid, centroid_no_neighbors);
}

void apply_rules(double &beta,
                 const std::vector<double> &c1,
                 const std::vector<double> &c2,
                 const std::vector<double> &current_position,
                 double dt,
                 double beta_min,
                 const double &betaD,
                 const std::vector<double> &goal,
                 double d1,
                 double &th,
                 double d2,
                 double d3,
                 double d4,
                 std::vector<double> &destinations,
                 std::vector<double> &c1_no_rotation)
{

    // Extract x, y components from current_position
    double current_j_x = current_position[0];
    double current_j_y = current_position[1];

    // first condition
    double dist_c1_c2 = sqrt(pow((c1[0] - c2[0]), 2) + pow((c1[1] - c2[1]), 2));
    if (dist_c1_c2 > d2 && sqrt(pow((current_j_x - c1[0]), 2) + pow((current_j_y - c1[1]), 2)) < d1)
    {
        beta = std::max(beta - dt, beta_min);
    }
    else
    {
        beta = beta - dt * (beta - betaD);
    }

    // second condition
    bool dist_c1_c2_d4 = dist_c1_c2 > d4;

    if (dist_c1_c2_d4 && sqrt(pow((current_j_x - c1[0]), 2) + pow((current_j_y - c1[1]), 2)) < d3)
    {
        th = std::min(th + dt, M_PI / 2);
    }
    else
    {
        th = std::max(0.0, th - dt);
    }

    // third condition
    if (th == M_PI / 2 && sqrt(pow((current_j_x - c1_no_rotation[0]), 2) + pow((current_j_y - c1_no_rotation[1]), 2)) > sqrt(pow((current_j_x - c1[0]), 2) + pow((current_j_y - c1[1]), 2)))
    {
        th = 0;
    }
    std::cout << (sqrt(pow((current_j_x - c1[0]), 2) + pow((current_j_y - c1[1]), 2)) < d3) << std::endl;
    // Compute the angle and new position
    double angle = atan2(goal[1] - current_j_y, goal[0] - current_j_x);
    double new_angle = angle - th;
    double distance = sqrt(pow((goal[0] - current_j_x), 2) + pow((goal[1] - current_j_y), 2));
    destinations[0] = current_j_x + distance * cos(new_angle);
    destinations[1] = current_j_y + distance * sin(new_angle);
}

int main()
{

    std::pair<double, double> robot_pos = {0.0, 0.0};
    double step_size = 0.075;
    double radius = 1.0;
    double beta = 0.2;
    std::vector<std::pair<double, double>> neighbors = {{0.5, -0.0}};
    std::pair<double, double> destination = {4.0, 0.0};

    std::vector<double> size_neighbors = {0.2, 0.2};
    double encumbrance = 0.2;
    double dt = 0.1;
    double beta_min = 0.07;
    double betaD = 0.2;
    const std::vector<double> goal = {destination.first, destination.second};
    double d1 = 0.1;
    double th;
    double d2 = 2 * encumbrance;
    double d3 = d1;
    double d4 = d2;
    // Call points_inside_circle function
    std::vector<std::pair<double, double>> points = points_inside_circle(robot_pos, radius, step_size);

    // Call find_closest_points function
    std::vector<std::pair<double, double>> closest_points = find_closest_points(robot_pos, points, neighbors);

    // Call account_encumbrance function
    std::vector<std::pair<double, double>> accounted_points = account_encumbrance(closest_points, robot_pos, neighbors, size_neighbors, encumbrance);

    // Call get_centroid function
    auto c1_c2 = get_centroid(robot_pos, radius, step_size, neighbors, size_neighbors, encumbrance, destination, beta);

    auto c1_c2_no_rotations = get_centroid(robot_pos, radius, step_size, neighbors, size_neighbors, encumbrance, {goal[0], goal[1]}, beta);

    std::vector<double> c1 = {c1_c2.first.first, c1_c2.first.second};
    std::vector<double> c2 = {c1_c2.second.first, c1_c2.second.second};

    std::vector<double> current_position(2); // Vector with two elements
    current_position[0] = robot_pos.first;   // Assigning first element
    current_position[1] = robot_pos.second;
    std::vector<double> destinations = {destination.first, destination.second};
    std::vector<double> c1_no_rotation = {c1_c2_no_rotations.first.first, c1_c2_no_rotations.first.second};

    // Call apply_rules function
    apply_rules(beta, c1, c2, current_position, dt, beta_min, betaD, goal, d1, th, d2, d3, d4, destinations, c1_no_rotation);

    // Output the updated values
    std::cout << "Updated beta: " << beta << std::endl;
    std::cout << "Updated destinations: ";
    for (auto value : destinations)
    {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    std::cout << "c1: ";
    for (const auto &value : c1)
    {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return 0;
}
