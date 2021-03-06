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
	num_particles = 50;
	default_random_engine gen;
	
	for(int i=0;i<num_particles;i++)
	{
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		double w = 1.0;
		p.weight = w;
		particles.push_back(p);
		weights.push_back(w);		
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	

	/*
	cout<<"\n\nPre-Prediction Check!!\n";
	for(int i=0;i<num_particles;i++)
	{
		cout<<"\nFor particle "<<i<<" the x is "<<particles[i].x;
	}
	*/
	
	// Calculate new x,y,theta
	double vel_by_yaw_rate = (velocity/yaw_rate);
	double vel_times_dt = velocity*delta_t;
	double yaw_rate_times_dt = yaw_rate*delta_t;
	default_random_engine gen;
	for(int i=0;i<num_particles;i++)
	{	
			
		double x_zero = particles[i].x;
		double y_zero = particles[i].y;
		double theta_zero = particles[i].theta;
		
		double x;
		double y;
		double theta;
		
		if((yaw_rate == 0)||(fabs(yaw_rate)<0.0001))
		{
			x = x_zero + vel_times_dt*cos(theta_zero);
			y = y_zero + vel_times_dt*sin(theta_zero);
			theta = theta_zero;
		}
		else 
		{
			x = x_zero + vel_by_yaw_rate*(sin(theta_zero + yaw_rate_times_dt) - sin(theta_zero));
			y = y_zero + vel_by_yaw_rate*(cos(theta_zero) - cos(theta_zero + yaw_rate_times_dt));
			theta = theta_zero + yaw_rate*delta_t;
		}
		
		// Initialize gaussians
		normal_distribution<double> dist_x(x, std_pos[0]);
		normal_distribution<double> dist_y(y, std_pos[1]);
		normal_distribution<double> dist_theta(theta, std_pos[2]);
		
		
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
	//cout<<"\n\n*~*~*~ Update Called *~*~*";
	for(int p=0; p<num_particles; p++)
	{
		//cout<<"\n\nStarting Particle: "<<p;
		double final_weight = 1.0;
		
		for(int o=0; o<observations.size(); o++)
		{
			//cout<<"\nParticle: "<<p<<" , Observation: "<<o;
					
			double transformed_x =   observations[o].x * cos(particles[p].theta)
								   - observations[o].y * sin(particles[p].theta)
								   + particles[p].x;
			
			double transformed_y =   observations[o].x * sin(particles[p].theta)
								   + observations[o].y * cos(particles[p].theta)
								   + particles[p].y;
								   
			
			double min_landmark_obs_distance = 200.0;
			int min_landmark_index = -1;
			
			for(int l=0; l<map_landmarks.landmark_list.size(); l++) {
				
				double landmark_particle_distance = dist(particles[p].x, particles[p].y,map_landmarks.landmark_list[l].x_f, map_landmarks.landmark_list[l].y_f);
														 
				if(landmark_particle_distance <= sensor_range)
				{ 
					double landmark_obs_distance = dist(transformed_x, transformed_y,
														map_landmarks.landmark_list[l].x_f, map_landmarks.landmark_list[l].y_f);
					
					if(landmark_obs_distance < min_landmark_obs_distance)
					{
						
						min_landmark_obs_distance = landmark_obs_distance;
						min_landmark_index = l;
						
					}
				}
				
			}
			 //cout<<"\nTransformed: ("<<transformed_x<<", "<<transformed_y<<")";
			 //cout<<"\nLandmark: ("<<map_landmarks.landmark_list[min_landmark_index].x_f<<", "<<map_landmarks.landmark_list[min_landmark_index].y_f<<")";
			if(min_landmark_index != -1)
			{
				
				double obs_probability = exp(-1*
												(
													(pow(transformed_x - map_landmarks.landmark_list[min_landmark_index].x_f, 2)
													/(2*pow(std_landmark[0],2)))
												  +
												    (pow(transformed_y - map_landmarks.landmark_list[min_landmark_index].y_f, 2)
													/(2*pow(std_landmark[1],2)))
												)
											)
											*1/(2 * M_PI * std_landmark[0] * std_landmark[1]);
				//cout<<"\nobs prob :"<<obs_probability;							
				final_weight *= obs_probability;
			}	
			
		}
		//cout<<"\nFinal Weight for Particle "<<p<<": "<<final_weight;
		particles[p].weight = final_weight;
		weights[p] = final_weight;
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
