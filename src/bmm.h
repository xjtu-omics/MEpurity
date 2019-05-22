#ifndef BMM_H
#define	BMM_H

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <sys/stat.h>

struct kmeans_output
{
	std::vector<float> _centers;
	std::vector<int> _clusters;
	float errorNum;
	kmeans_output() {}
	kmeans_output(const kmeans_output &output) {
		_centers = output._centers;
		_clusters = output._clusters;
		errorNum = output.errorNum;
	}
};

struct bmm_outputs
{
	std::vector<std::vector<float> > ln_rho;
	std::vector<float> E_lnu;
	std::vector<float> E_lnv;
	std::vector<float> E_lnpi;
	std::vector<float> E_quadratic_u;
	std::vector<float> E_quadratic_v;
	std::vector<float> u_bar;
	std::vector<float> v_bar;
	int iteration_time;
	bmm_outputs() {}
	bmm_outputs(const bmm_outputs &output0) {
		ln_rho = output0.ln_rho;
		E_lnu = output0.E_lnu;
		E_lnpi = output0.E_lnpi;
		E_quadratic_u = output0.E_quadratic_u;
		E_quadratic_v = output0.E_quadratic_v;
		u_bar = output0.u_bar;
		v_bar = output0.v_bar;
		iteration_time = output0.iteration_time;
	}
};

struct bmm_results
{
	std::map<float, int> cluster_result;
	std::vector<float> cluster_mean;
	std::vector<int> cluster_num;
	bmm_results() {}
	~bmm_results() {}
	bmm_results(const bmm_results &bmm_result) {
		cluster_result = bmm_result.cluster_result;
		cluster_mean = bmm_result.cluster_mean;
		cluster_num = bmm_result.cluster_num;
	}
};

struct hyperparameter {
	std::vector<float> mu0;
	std::vector<float> alpha0;
	std::vector<float> nu0;
	std::vector<float> beta0;
	std::vector<float> c0;
	hyperparameter() {}
	~hyperparameter() {}
	hyperparameter(const hyperparameter &hyperparameters) {
		mu0 = hyperparameters.mu0;
		alpha0 = hyperparameters.alpha0;
		nu0 = hyperparameters.nu0;
		beta0 = hyperparameters.beta0;
		c0 = hyperparameters.c0;
	}

};

struct parameter
{
	std::vector<float> mu1;
	std::vector<float> alpha1;
	std::vector<float> nu1;
	std::vector<float> beta1;
	std::vector<float> c1;
	std::vector<float> epi1;
	std::vector<std::vector<float> > r1;
	std::vector<float> kmeans_centers;
	std::vector<int> kmeans_clusters;
	parameter() {}
	~parameter() {}
	parameter(const parameter &parameters){
		mu1 = parameters.mu1;
		alpha1 = parameters.mu1;
		nu1 = parameters.nu1;
		beta1 = parameters.beta1;
		c1 = parameters.c1;
		epi1 = parameters.epi1;
		r1 = parameters.r1;
		kmeans_centers = parameters.kmeans_centers;
		kmeans_clusters = parameters.kmeans_clusters;
	}
};

std::ostream & operator<<(std::ostream &out, std::vector<float> &obj);
std::ofstream & operator<<(std::ofstream &out, bmm_results &obj);
void init_bmm_parameters(parameter& parameters, hyperparameter hyperparameters, std::vector<float> sample, int class_number) ;
void init_bmm_parameters_text(parameter& parameters, hyperparameter hyperparameters, kmeans_output kmeans_outputs, std::vector<float> sample, int class_number) ;
void init_bmm_hyperparameters(hyperparameter& hyperparameters, int class_count);
bmm_outputs bmm_one_run(std::vector<float> sample, int class_number, parameter& parameters, hyperparameter &hyperparameters , float threshold = 0.0001 , unsigned int max_iter = 10000) ;
bmm_results bmm_function(std::vector<float> sample, int class_number, parameter& parameters, hyperparameter &hyperparameters, float cut_threshold = 0.01, unsigned int max_iter = 10000, float threshold = 0.0001) ;
void show_kmeans_output(kmeans_output kmoutput) ;
std::vector<float> read_file(const char* filename);
std::vector<float> read_file2(const char* filename);
void load_kmeans_parameters(const char* filename, kmeans_output &kmeans_outputs);
kmeans_output my_kmeans_little(std::vector<float> sample);
#endif