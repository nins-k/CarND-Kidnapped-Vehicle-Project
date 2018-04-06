/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	// Initialize gaussians
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	// Initialize particles
	num_particles = 5;
	default_random_engine gen;
	
	for(int i=0;i<num_particles;i++)
	{
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		particles.push_back(p);
		double w = 1;
		p.weight = w;
		weights.push_back(w);
	}
	
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	double x;
	double y;
	double theta;
	
	// Calculate new x,y,theta
	for(int i=0;i<num_particles;i++)
	{	
		double theta_zero = particles[i].theta;
		double x_zero = particles[i].x;
		double y_zero = particles[i].y;
		
		if(yaw_rate == 0)
		{
			x = x_zero + velocity*delta_t*cos(theta_zero);
			y = y_zero + velocity*delta_t*sin(theta_zero);
			theta = theta_zero;
		}
		else 
		{
			x = x_zero + (velocity/yaw_rate)*(sin(theta_zero + yaw_rate*delta_t) - sin(theta_zero));
			y = y_zero + (velocity/yaw_rate)*(cos(theta_zero) - cos(theta_zero + yaw_rate*delta_t));
			theta = theta_zero + yaw_rate*delta_t;
		}
		
		// Initialize gaussians
		normal_distribution<double> dist_x(x, std_pos[0]);
		normal_distribution<double> dist_y(y, std_pos[1]);
		normal_distribution<double> dist_theta(theta, std_pos[2]);
		default_random_engine gen;
		
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}
	
	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	// For each particle ~
	for(int t=0;t<num_particles;t++)
	{
		Particle particle = particles[t];
		
		// Obtain a vector of observations in map co-ordinates
		std::vector<LandmarkObs> map_coord_observations;
		for(int i=0;i<observations.size();i++)
		{
			LandmarkObs car_obs = observations[i];
			LandmarkObs trans_obs;
			
			trans_obs.x = car_obs.x*cos(particle.theta) - car_obs.y*sin(particle.theta) + particle.x;
			trans_obs.y = car_obs.x*sin(particle.theta) + car_obs.y*cos(particle.theta) + particle.y;
			trans_obs.observed_distance = dist(particle.x, particle.y, trans_obs.y, trans_obs.y);
			map_coord_observations.push_back(trans_obs);
		}
		
		// Obtain a vector of landmarks that are within the sensor range
		std::vector<int> inrange_landmarks;
		for(int i=0;i<map_landmarks.landmark_list.size();i++)
		{
			Map::single_landmark_s landmark = map_landmarks.landmark_list[i];
			float expected_distance = dist(particle.x, particle.y, landmark.x_f, landmark.y_f);
			if(sensor_range<=expected_distance)
			{
				inrange_landmarks.push_back(i);
			}
		}
		
		// Obtain index of associated landmarks for each observation
		std::vector<int> associated_landmarks;	
		for(int i=0; i<map_coord_observations.size(); i++)
		{
			double min_dist = 1000000.0;
			double associated_landmark_index = -1;
			
			LandmarkObs obs = map_coord_observations[i];
			
			for(int j=0; j<inrange_landmarks.size(); j++)
			{
				Map::single_landmark_s landmark = map_landmarks.landmark_list[inrange_landmarks[j]];
				float expected_distance = dist(obs.x, obs.y, landmark.x_f, landmark.y_f);
				if(expected_distance<min_dist)
				{
					associated_landmark_index = inrange_landmarks[j];
					min_dist = expected_distance;
				}
			}
			cout<<"\n\nmin dist for particle "<<t<< " is "<<min_dist;
			cout<<"\nparticle x is "<<particle.x<<", particle y is "<<particle.y;
			
			LandmarkObs lm = map_coord_observations[associated_landmark_index];
			cout<<"\nobs x is "<<obs.x<<", obs y is "<<obs.y;
			cout<<"\nAL x is " <<lm.x<<", AL y is "<<lm.y;
			associated_landmarks.push_back(associated_landmark_index);
		}
		
		// Calculate the particle's final weight
		double final_weight = 1.0;
		
		for(int k=0; k<map_coord_observations.size(); k++)
		{
			LandmarkObs obs = map_coord_observations[k];
			Map::single_landmark_s landmark = map_landmarks.landmark_list[associated_landmarks[k]];
			
			double denominator = 1.0/(2.0*M_PI*std_landmark[0]*std_landmark[1]);
			double exp_term_x = (obs.x - landmark.x_f)*(obs.x - landmark.x_f)/(2.0*std_landmark[0]*std_landmark[0]);
			double exp_term_y = (obs.y - landmark.y_f)*(obs.y - landmark.y_f)/(2.0*std_landmark[1]*std_landmark[1]);
			//std::cout<<"exp_term_x: "<<exp_term_x<<endl;
			//std::cout<<"exp_term_y: "<<exp_term_y<<endl;
			//std::cout<<"Final weight pre exp: "<<final_weight<<endl;
			final_weight = final_weight* (denominator * exp(- exp_term_x - exp_term_y));
			//std::cout<<"Final weight post exp: "<<final_weight<<endl;
		}
		particle.weight = final_weight;
		weights[t] = final_weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());
	
	std::vector<Particle> resample_particles;
	
	for(int i=0; i<num_particles; i++)
	{
		resample_particles.push_back(particles[distribution(gen)]);
	}
	
	particles = resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
