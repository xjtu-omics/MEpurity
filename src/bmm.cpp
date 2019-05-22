#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <sys/stat.h>
#include <utility>
#include <numeric>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include "bmm.h"
#include "kmeans.h"
using namespace std;
using namespace boost;


std::ostream & operator<<(std::ostream &out, std::vector<float> &obj) {
	for (int i = 0; i < obj.size(); i++) {
		out << obj[i] << '\t';
	}
	out << endl;
	return out;
}

std::ofstream & operator<<(std::ofstream &out, bmm_results &obj) {
	std::map<float, int> sort_results;
	for(int i = 0; i < obj.cluster_mean.size(); i++){
		sort_results.insert(make_pair(obj.cluster_mean[i], obj.cluster_num[i]));
	}
	for(std::map<float, int>::iterator it = sort_results.begin(); it != sort_results.end(); it++){
		out << '\t' << it->first << '\t' << it->second;
	}
	out << endl;
	return out;
}

kmeans_output my_kmeans1(std::vector<float> data, std::vector<int> lables, int k_class = 10, unsigned int max_iter = 100) {
	kmeans_output outputtt;
	double *a;
	double *c;
	int i;
	int *ic1;
	int ifault;
	std::ifstream input;
	int iter;
	int j;
	int k = k_class;
	int m = data.size();
	int n = 1;
	int *nc;
	int nc_sum;
	double *wss;
	double wss_sum;
	a = new double[m*n];
	c = new double[k*n];
	ic1 = new int[m];
	nc = new int[k];
	wss = new double[k];
	for (i = 0; i < m; i++) {
		a[i] = data[i];
	}
	input.close();
	for (i = 1; i <= k; i++){
		c[i - 1] = a[lables[i-1]];
	}
	iter = max_iter;
	kmns(a, m, n, c, k, ic1, nc, iter, wss, &ifault);
	if (ifault != 0 && ifault != 3)
	{
		outputtt.errorNum = 10000;
		delete[] a;
		delete[] c;
		delete[] ic1;
		delete[] nc;
		delete[] wss;
		return outputtt;
	}
	nc_sum = 0;
	wss_sum = 0.0;
	for (i = 1; i <= k; i++)
	{
		wss_sum = wss_sum + wss[i - 1];
	}
	vector<float> centers(k, 0);
	vector<int> clusters(m, 0);
	for (i = 0; i < k; i++) {
		centers[i] = c[i];
	}
	for (i = 0; i < m; i++) {
		clusters[i] = ic1[i]-1;
	}
	delete[] a;
	delete[] c;
	delete[] ic1;
	delete[] nc;
	delete[] wss;
	outputtt._centers = centers;
	outputtt._clusters = clusters;
	outputtt.errorNum = wss_sum;
	//cout<<wss_sum<<endl;
	return outputtt;
}

kmeans_output my_kmeans(std::vector<float> sample, int k_class=10, unsigned int nstart=1000, unsigned int max_iter=100) {
	kmeans_output ansmeans;
	ansmeans.errorNum = 10000;
	srand((unsigned)time(NULL));
	int sizez=sample.size();
	for (int i = 0; i < nstart; i++) {
		vector<int> lables;
		for(int j = 0; j < k_class; j++){
			lables.push_back(rand()% sizez);
		}
		kmeans_output thisansmeans = my_kmeans1(sample,lables, k_class, max_iter);
		if (thisansmeans._clusters.size() > 0 && thisansmeans.errorNum < ansmeans.errorNum) {
			ansmeans = thisansmeans;
		}
	}
	return ansmeans;
}

kmeans_output my_kmeans_little(std::vector<float> sample){
	float min_beta = sample[sample.size() - 1];
	float max_beta = sample[0];
	kmeans_output ansmeans;
	std::vector<float> centers;
	std::vector<int> clusters;
	if(sample.size() == 2){
		centers.push_back(sample[0]);
		centers.push_back(sample[1]);
		clusters.push_back(0);
		clusters.push_back(1);
	}
	else{
		if(sample[1] - sample[0] < sample[2] - sample[1]){
			centers.push_back((sample[0] + sample[1]) / 2);
			centers.push_back(sample[2]);
			clusters.push_back(0);
			clusters.push_back(0);
			clusters.push_back(1);
		}
		else{
			centers.push_back(sample[0]);
			centers.push_back((sample[1] + sample[2]) / 2);
			clusters.push_back(0);
			clusters.push_back(1);
			clusters.push_back(1);
		}
	}
	ansmeans._centers = centers;
	ansmeans._clusters = clusters;
	return ansmeans;
}

void init_bmm_parameters(parameter& parameters, hyperparameter hyperparameters, std::vector<float> sample, int class_number) {
	kmeans_output origion_sites = my_kmeans(sample,class_number);
	if(origion_sites._centers.size() == 0){
		cerr << "Fail to initial BMM parameters. Please reduce the number of initial clusters.\n";
		exit(1);
	}
	parameters.kmeans_centers = origion_sites._centers;
	parameters.kmeans_clusters = origion_sites._clusters;
	int data_size = sample.size();
	std::vector<float> r_0(class_number, 0);
	std::vector<std::vector<float> >r_matrix;
	for (int i = 0; i < data_size; i++) {
		std::vector<float> r_1(class_number,0);
		r_1[parameters.kmeans_clusters[i]] = 1;
		r_matrix.push_back(r_1);
	}
	parameters.r1=r_matrix;
	std::vector<float> u_bar1 = r_0;
	std::vector<float> v_bar1 = r_0;
	std::vector<float> e_lnu1 = r_0;
	std::vector<float> e_lnv1 = r_0;
	std::vector<float> e_lnpi1 = r_0;
	std::vector<float> e_pi1 = r_0;
	std::vector<float> E_quadratic_u1 = r_0;
	std::vector<float> E_quadratic_v1 = r_0;
	std::vector<float> lg_u_v1 = r_0;
	std::vector<float> lg_u1 = r_0;
	std::vector<float> lg_v1 = r_0;
	std::vector<float> dig_u_v1 = r_0;
	std::vector<float> dig_u1 = r_0;
	std::vector<float> dig_v1 = r_0;
	std::vector<float> trig_u_v1 = r_0;
	std::vector<float> trig_u1 = r_0;
	std::vector<float> trig_v1 = r_0;
	std::vector<float> E_lnv_logvbar1 = r_0;
	std::vector<float> E_lnu_logvbar1 = r_0;
	std::vector<float> v_trig_E_log1 = r_0;
	std::vector<float> u_trig_E_log1 = r_0;
	std::vector<float> r_colsums1 = r_0;
	parameters.mu1 = hyperparameters.mu0;
	parameters.nu1 = hyperparameters.nu0;
	parameters.beta1 = hyperparameters.beta0;
	parameters.alpha1 = hyperparameters.alpha0;
	parameters.c1 = r_0;
	parameters.epi1 = r_0;
	float c_sum = 0;
	for (int i = 0; i < class_number; i++) {
		parameters.c1[i] = hyperparameters.c0[i];
		for (int j = 0; j < data_size; ++j) {
			parameters.c1[i] += parameters.r1[j][i];
		}
		c_sum += parameters.c1[i];
	}
	for (int i = 0; i < class_number; i++) {
		parameters.alpha1[i] = parameters.mu1[i] * parameters.beta1[i] * (1 - parameters.kmeans_centers[i]) / parameters.kmeans_centers[i] * parameters.nu1[i];
		u_bar1[i] = parameters.mu1[i] / parameters.alpha1[i];
		v_bar1[i] = parameters.nu1[i] / parameters.beta1[i];
		e_lnu1[i] = boost::math::digamma(parameters.mu1[i]) - log(parameters.alpha1[i]);
		e_lnv1[i] = boost::math::digamma(parameters.nu1[i]) - log(parameters.beta1[i]);
		e_lnpi1[i] = boost::math::digamma(parameters.c1[i]) - boost::math::digamma(c_sum);
		E_quadratic_u1[i] = (boost::math::digamma(parameters.mu1[i]) - log(parameters.mu1[i]))*(boost::math::digamma(parameters.mu1[i]) - log(parameters.mu1[i])) + boost::math::trigamma(parameters.mu1[i]);
		E_quadratic_v1[i] = (boost::math::digamma(parameters.nu1[i]) - log(parameters.nu1[i]))*(boost::math::digamma(parameters.nu1[i]) - log(parameters.nu1[i])) + boost::math::trigamma(parameters.nu1[i]);
		lg_u_v1[i] = boost::math::lgamma(u_bar1[i] + v_bar1[i]);
		lg_u1[i] = boost::math::lgamma(u_bar1[i]);
		lg_v1[i] = boost::math::lgamma(v_bar1[i]);
		dig_u_v1[i] = boost::math::digamma(u_bar1[i] + v_bar1[i]);
		dig_u1[i] = boost::math::digamma(u_bar1[i]);
		dig_v1[i] = boost::math::digamma(v_bar1[i]);
		trig_u_v1[i] = boost::math::trigamma(u_bar1[i] + v_bar1[i]);
		trig_u1[i] = boost::math::trigamma(u_bar1[i]);
		trig_v1[i] = boost::math::trigamma(v_bar1[i]);
		E_lnu_logvbar1[i] = e_lnu1[i] - log(u_bar1[i]);
		E_lnv_logvbar1[i] = e_lnv1[i] - log(v_bar1[i]);
		parameters.alpha1[i] = hyperparameters.alpha0[i];
		for (int j = 0; j < data_size; ++j) {
			parameters.alpha1[i] -= parameters.r1[j][i] * log(sample[j]);
			parameters.beta1[i] -= parameters.r1[j][i] * log(1 - sample[j]);
			r_colsums1[i] += parameters.r1[j][i];
		}
		v_trig_E_log1[i] = v_bar1[i] * trig_u_v1[i] * E_lnv_logvbar1[i];
		u_trig_E_log1[i] = u_bar1[i] * trig_u_v1[i] * E_lnu_logvbar1[i];
		parameters.mu1[i] += r_colsums1[i] * u_bar1[i] * (dig_u_v1[i] - dig_u1[i] + v_trig_E_log1[i]);
		parameters.nu1[i] += r_colsums1[i] * v_bar1[i] * (dig_u_v1[i] - dig_v1[i] + u_trig_E_log1[i]);
	}
	std::cout << "Parameters calculate over." << endl;
}


void init_bmm_hyperparameters(hyperparameter& hyperparameters, int class_count) {
	std::vector<float> mu0(class_count, 1);
	std::vector<float> nu0(class_count, 1);
	std::vector<float> alpha0(class_count, 0.005);
	std::vector<float> beta0(class_count, 0.005);
	std::vector<float> c0(class_count, 0.001);
	hyperparameters.mu0 = mu0;
	hyperparameters.nu0 = nu0;
	hyperparameters.alpha0 = alpha0;
	hyperparameters.beta0 = beta0;
	hyperparameters.c0 = c0;
	std::cout << "Hyperparameters load over." << endl;
}



bmm_outputs bmm_one_run(std::vector<float> sample, int class_number, parameter& parameters, hyperparameter &hyperparameters, float threshold, unsigned int max_iter) {
	std::vector<float> r_0(class_number, 0);
	std::vector<float> u_bar = r_0;
	std::vector<float> v_bar = r_0;
	std::vector<float> E_pi_prev = parameters.epi1;
	std::vector<float> e_lnu = r_0;
	std::vector<float> e_lnv = r_0;
	std::vector<float> e_lnpi = r_0;
	std::vector<float> E_quadratic_u = r_0;
	std::vector<float> E_quadratic_v = r_0;
	std::vector<float> lg_u_v = r_0;
	std::vector<float> lg_u = r_0;
	std::vector<float> lg_v = r_0;
	std::vector<float> dig_u_v = r_0;
	std::vector<float> dig_u = r_0;
	std::vector<float> dig_v = r_0;
	std::vector<float> trig_u_v = r_0;
	std::vector<float> trig_u = r_0;
	std::vector<float> trig_v = r_0;
	std::vector<float> E_lnv_logvbar = r_0;
	std::vector<float> E_lnu_logubar = r_0;
	std::vector<float> u_v_tmp = r_0;
	std::vector<float> v_trig_E_log = r_0;
	std::vector<float> u_trig_E_log = r_0;
	std::vector<std::vector<float> > u_logx;
	std::vector<std::vector<float> > v_logx;

	std::vector<std::vector<float> > r_matrix;

	bmm_outputs bmm_output;


	int data_size = sample.size();
	int iteration_time = 0;
	float c_sum = 0;
	float sum_c0 = 0;
	std::vector<float> r_00(data_size, 0);
	for (int i = 0; i < class_number; i++) {
		sum_c0 += hyperparameters.c0[i];
		u_bar[i] = parameters.mu1[i] / parameters.alpha1[i];
		v_bar[i] = parameters.nu1[i] / parameters.beta1[i];
		c_sum += parameters.c1[i];
		u_logx.push_back(r_00);
		v_logx.push_back(r_00);
	}
	for (int i = 0; i < data_size; i++) {
		r_matrix.push_back(r_0);
	}

	while (1) {
		std::vector<std::vector<float> > ln_rho;
		for (int i = 0; i < class_number; i++) {
			ln_rho.push_back(r_00);
		}
		std::vector<float> one_vector_sum = r_00;
		std::vector<float> mutliple_result1 = r_0;
		std::vector<float> mutliple_result2 = r_0;
		std::vector<float> cow_sum = r_0;
		for (int i = 0; i < class_number; i++) {
			e_lnu[i] = boost::math::digamma(parameters.mu1[i]) - log(parameters.alpha1[i]);
			e_lnv[i] = boost::math::digamma(parameters.nu1[i]) - log(parameters.beta1[i]);
			e_lnpi[i] = boost::math::digamma(parameters.c1[i]) - boost::math::digamma(c_sum);
			E_quadratic_u[i] = (boost::math::digamma(parameters.mu1[i]) - log(parameters.mu1[i]))*(boost::math::digamma(parameters.mu1[i]) - log(parameters.mu1[i])) + boost::math::trigamma(parameters.mu1[i]);
			E_quadratic_v[i] = (boost::math::digamma(parameters.nu1[i]) - log(parameters.nu1[i]))*(boost::math::digamma(parameters.nu1[i]) - log(parameters.nu1[i])) + boost::math::trigamma(parameters.nu1[i]);
			lg_u_v[i] = boost::math::lgamma(u_bar[i] + v_bar[i]);
			lg_u[i] = boost::math::lgamma(u_bar[i]);
			lg_v[i] = boost::math::lgamma(v_bar[i]);
			dig_u_v[i] = boost::math::digamma(u_bar[i] + v_bar[i]);
			dig_u[i] = boost::math::digamma(u_bar[i]);
			dig_v[i] = boost::math::digamma(v_bar[i]);
			trig_u_v[i] = boost::math::trigamma(u_bar[i] + v_bar[i]);
			trig_u[i] = boost::math::trigamma(u_bar[i]);
			trig_v[i] = boost::math::trigamma(v_bar[i]);
			E_lnu_logubar[i] = e_lnu[i] - log(u_bar[i]);
			E_lnv_logvbar[i] = e_lnv[i] - log(v_bar[i]);
			u_v_tmp[i] = lg_u_v[i] - lg_u[i] - lg_v[i] + u_bar[i] * (dig_u_v[i] - dig_u[i])*E_lnu_logubar[i] + v_bar[i] * (dig_u_v[i] - dig_v[i])*E_lnv_logvbar[i] + 0.5*u_bar[i] * u_bar[i] * (trig_u_v[i] - trig_u[i])\
				*E_quadratic_u[i] + 0.5*v_bar[i] * v_bar[i] * (trig_u_v[i] - trig_v[i])*E_quadratic_v[i] + u_bar[i] * v_bar[i] * trig_u_v[i] * E_lnu_logubar[i] * E_lnv_logvbar[i];
			float u_bar_1 = u_bar[i] - 1;
			float v_bar_1 = v_bar[i] - 1;
			for (int j = 0; j < data_size; j++) {
				u_logx[i][j] = log(sample[j])*u_bar_1;
				v_logx[i][j] = log(1 - sample[j])*v_bar_1;
				ln_rho[i][j] = e_lnpi[i] + u_v_tmp[i] + u_logx[i][j] + v_logx[i][j];//Is anything wrong?				
			}
		}
		for (int j = 0; j < data_size; j++) {
			float row_sum = 0;
			float max_rho = -10000000;
			for (int i = 0; i < class_number; i++) {
				if (ln_rho[i][j] > max_rho) {
					max_rho = ln_rho[i][j];
				}
			}
			for (int i = 0; i < class_number; i++) {
				one_vector_sum[j] += std::exp(ln_rho[i][j] - max_rho);

			}
			for (int i = 0; i < class_number; i++) {
				r_matrix[j][i] = std::exp(ln_rho[i][j] - log(one_vector_sum[j]) - max_rho) + 1e-9;//是否要转为double？
				cow_sum[i] += r_matrix[j][i];
				mutliple_result1[i] += r_matrix[j][i] * log(sample[j]);
				mutliple_result2[i] += r_matrix[j][i] * log(1 - sample[j]);
			}
		}
		for (int i = 0; i < class_number; i++) {
			parameters.epi1[i] = (hyperparameters.c0[i] + cow_sum[i]) / (sum_c0 + data_size);
			parameters.alpha1[i] = hyperparameters.alpha0[i] - mutliple_result1[i];
			parameters.beta1[i] = hyperparameters.beta0[i] - mutliple_result2[i];
			v_trig_E_log[i] = v_bar[i] * trig_u_v[i] * E_lnv_logvbar[i];
			parameters.mu1[i] = hyperparameters.mu0[i] + cow_sum[i] * u_bar[i] * (dig_u_v[i] - dig_u[i] + v_trig_E_log[i]);
			u_trig_E_log[i] = u_bar[i] * trig_u_v[i] * E_lnu_logubar[i];
			parameters.nu1[i] = hyperparameters.nu0[i] + cow_sum[i] * v_bar[i] * (dig_u_v[i] - dig_v[i] + u_trig_E_log[i]);
			parameters.c1[i] = hyperparameters.c0[i] + cow_sum[i];
			u_bar[i] = parameters.mu1[i] / parameters.alpha1[i];
			v_bar[i] = parameters.nu1[i] / parameters.beta1[i];
		}
		int more_count = 0;
		for (int i = 0; i < class_number; i++) {
			if (std::abs(parameters.epi1[i] - E_pi_prev[i]) < threshold) {
				more_count++;
			}
		}
		if (more_count >= class_number) {
			bmm_output.E_lnpi = e_lnpi;
			bmm_output.E_lnu = e_lnu;
			bmm_output.E_lnv = e_lnv;
			bmm_output.E_quadratic_u = E_quadratic_u;
			bmm_output.E_quadratic_v = E_quadratic_v;
			bmm_output.ln_rho = ln_rho;
			bmm_output.u_bar = u_bar;
			bmm_output.v_bar = v_bar;
			bmm_output.iteration_time = iteration_time;
			return bmm_output;
		}
		E_pi_prev = parameters.epi1;
		++iteration_time;
		if (iteration_time >= max_iter) {
			bmm_output.E_lnpi = e_lnpi;
			bmm_output.E_lnu = e_lnu;
			bmm_output.E_lnv = e_lnv;
			bmm_output.E_quadratic_u = E_quadratic_u;
			bmm_output.E_quadratic_v = E_quadratic_v;
			bmm_output.ln_rho = ln_rho;
			bmm_output.u_bar = u_bar;
			bmm_output.v_bar = v_bar;
			bmm_output.iteration_time = iteration_time;
			return bmm_output;
		}
	}
}



bmm_results bmm_function(std::vector<float> sample, int class_number, parameter& parameters, hyperparameter &hyperparameters, float cut_threshold, unsigned int max_iter, float threshold) {
	std::cout<<"Bmm model begin."<<endl;
	int total_iteration_count = 0;
	int data_size = sample.size();
	for (int i = 0; i < class_number; i++) {
		parameters.epi1[i] = 1 / class_number;
	}
	std::vector<float> r_0(class_number, 0);
	std::vector<std::vector<float> > ln_rho_0;
	std::vector<float> mu_0, alpha_0, nu_0, beta_0, c_0, epi_0, e_lnu_0, e_lnv_0, e_lnpi_0, E_quadratic_u_0, E_quadratic_v_0, u_bar_0, v_bar_0;
	bool do_cut_sites = 0;
	int cluster_number = class_number;
	int preclusternum = cluster_number;
	parameter preparameters;
	int least_class_num;
	if(class_number>10){
		least_class_num = 0;
	}
	else if(data_size > 500){
		least_class_num = 10;
	}
	else{
		least_class_num = class_number/2;
	}
	while (1) {
		bmm_outputs bmm_output=bmm_one_run(sample, cluster_number, parameters, hyperparameters, threshold);
		mu_0 = parameters.mu1;
		alpha_0 = parameters.alpha1;
		nu_0 = parameters.nu1;
		beta_0 = parameters.beta1;
		c_0 = parameters.c1;
		epi_0 = parameters.epi1;
		e_lnu_0 = bmm_output.E_lnu;
		e_lnv_0 = bmm_output.E_lnv;
		e_lnpi_0 = bmm_output.E_lnpi;
		ln_rho_0 = bmm_output.ln_rho;
		E_quadratic_u_0 = bmm_output.E_quadratic_u;
		E_quadratic_v_0 = bmm_output.E_quadratic_v;
		u_bar_0 = bmm_output.u_bar;
		v_bar_0 = bmm_output.v_bar;
		vector<int> vector_lables;
		if (cluster_number >= least_class_num) {
			for (int i = 0; i < cluster_number; i++) {
				if (epi_0[i] > cut_threshold) {
					vector_lables.push_back(i);
				}
			}
		}
		if(cluster_number < least_class_num){
			parameters = preparameters;
			cluster_number = preclusternum;
			break;
		}
		if(vector_lables.size()==0){
			break;
		}
		preparameters = parameters;
		preclusternum = cluster_number;
		if (vector_lables.size() < cluster_number) {
			std::vector<float> epi_1, e_lnpi_1, c_1, c_0_1, mu_1, nu_1, mu_0_1, nu_0_1, alpha_1, beta_1, alpha_0_1, beta_0_1;
			std::vector<std::vector<float> > r_1, ln_rho_1;
			for (int i = 0; i < vector_lables.size(); i++) {
				cluster_number = vector_lables.size();
				int at_location = vector_lables[i];
				epi_1.push_back(epi_0[at_location]);
				e_lnpi_1.push_back(e_lnpi_0[at_location]);				
				c_1.push_back(c_0[at_location]);
				c_0_1.push_back(hyperparameters.c0[at_location]);
				ln_rho_1.push_back(ln_rho_0[at_location]);
				mu_1.push_back(parameters.mu1[at_location]);//从这开始改
				nu_1.push_back(parameters.nu1[at_location]);
				mu_0_1.push_back(hyperparameters.mu0[at_location]);
				nu_0_1.push_back(hyperparameters.nu0[at_location]);
				alpha_1.push_back(parameters.alpha1[at_location]);
				beta_1.push_back(parameters.beta1[at_location]);
				alpha_0_1.push_back(hyperparameters.alpha0[at_location]);
				beta_0_1.push_back(hyperparameters.beta0[at_location]);
			}
			for (int j = 0; j < data_size; j++) {
				std::vector<float> one_line;
				for (int i = 0; i < vector_lables.size(); i++) {
					one_line.push_back(parameters.r1[j][vector_lables[i]]);
				}
				r_1.push_back(one_line);
			}
			std::vector<float> one_vector_sum(data_size, 0);

			for (int j = 0; j < data_size; j++) {
				float row_sum = 0;
				float max_rho = -10000000;
				for (int i = 0; i < cluster_number; i++) {
					if (ln_rho_1[i][j] > max_rho) {
						max_rho = ln_rho_1[i][j];
					}
				}
				for (int i = 0; i < cluster_number; i++) {
					one_vector_sum[j] += std::exp(ln_rho_1[i][j] - max_rho);

				}
				for (int i = 0; i < cluster_number; i++) {
					r_1[j][i] = std::exp(ln_rho_1[i][j] - log(one_vector_sum[j]) - max_rho);//是否要转为double？
				}
			}
			hyperparameters.alpha0 = alpha_0_1;
			hyperparameters.beta0 = beta_0_1;
			hyperparameters.c0 = c_0_1;
			hyperparameters.mu0 = mu_0_1;
			hyperparameters.nu0 = nu_0_1;
			parameters.alpha1 = alpha_1;
			parameters.beta1 = beta_1;
			parameters.c1 = c_1;
			parameters.epi1 = epi_1;
			parameters.mu1 = mu_1;
			parameters.nu1 = nu_1;
			parameters.r1 = r_1;
		}
		else {
			break;
		}
	}
	std::vector<float> cluster_mean;
	std::map<float, int> cluster_result;
	std::vector<int> cluster_sites_num(cluster_number, 0);
	std::vector<float> cluster_sites_sum(cluster_number, 0);
	float max_c_value = -1;
	int class_type = -1;
	for (int i = 0; i < data_size; i++) {
		max_c_value = -1;
		class_type = -1;
		for (int j = 0; j < cluster_number; j++) {
			if (parameters.r1[i][j] > max_c_value) {
				max_c_value = parameters.r1[i][j];
				class_type = j;
			}
		}
		cluster_sites_num[class_type]++;
		cluster_sites_sum[class_type] += sample[i];
		cluster_result.insert(make_pair(sample[i], class_type));
	}
	for (int i = 0; i < cluster_number; i++) {
		if( cluster_sites_num[i]!=0)
		cluster_mean.push_back(cluster_sites_sum[i] / cluster_sites_num[i]);
	}
	bmm_results bmm_result;
	bmm_result.cluster_mean = cluster_mean;
	bmm_result.cluster_result = cluster_result;
	bmm_result.cluster_num = cluster_sites_num;
	return bmm_result;
}



std::vector<float> read_file(const char* filename) {
	std::vector<float> sites_beta;
	std::filebuf *pbuf;
	std::ifstream filestr;
	long size;
	char *buffer;
	filestr.open(filename, ios::binary);
	pbuf = filestr.rdbuf();
	size = pbuf->pubseekoff(0, ios::end, ios::in);
	pbuf->pubseekpos(0, ios::in);
	buffer = new char[size];
	pbuf->sgetn(buffer, size);

	filestr.close();
	char *p = buffer;
	char* q = p;
	float beta_value;
	while ((p - buffer) < size) {
		while (*p != '\t' && (p - buffer)<size) { p++; }
		p++;
		q = p;
		while (*p != '\n' && (p - buffer)<size) { p++; }
		std::string num;
		while (q < p) {
			num += *q;
			q++;
		}
		beta_value = atof(num.c_str());
		p++;
		if(beta_value>0&&beta_value<1){
			sites_beta.push_back(beta_value);
		}
	}
	delete[] buffer;
	std::cout << "File load over." << endl;
	return sites_beta;
}

std::vector<float> read_file2(const char* filename) {
	std::vector<float> sites_beta;
	std::filebuf *pbuf;
	std::ifstream filestr;
	long size;
	char *buffer;
	filestr.open(filename, ios::binary);
	pbuf = filestr.rdbuf();
	size = pbuf->pubseekoff(0, ios::end, ios::in);
	pbuf->pubseekpos(0, ios::in); 
	buffer = new char[size];
	pbuf->sgetn(buffer, size);
	filestr.close();
	char *p = buffer;
	char* q = p;
	float beta_value;
	while ((p - buffer) < size) {
		while (*p != '\t' && (p - buffer)<size) { p++; }
		p++;
		q = p;
		while (*p != '\t' && (p - buffer)<size) { p++; }
		std::string num;
		while (q < p) {
			num += *q;
			q++;
		}
		beta_value = atof(num.c_str());
		p++;
		if(beta_value>0&&beta_value<1){
			sites_beta.push_back(beta_value);
		}
		while (*p != '\n' && (p - buffer)<size) { p++; }
		p++;
	}
	delete[] buffer;
	std::cout << "File load over." << endl;
	return sites_beta;
}




