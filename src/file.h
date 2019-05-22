#ifndef FILE_H
#define FILE_H

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
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

void run_all_files(const char* filepath, std::map<std::string, std::vector<float> >& distances);
void sort_and_output_map(const char* filename, std::map<std::string, float>& aftermap, float upper = 0);
void calculate_variance(std::map<std::string, std::vector<float> > distances, std::map<std::string, std::vector<float> >& viriation);
void calculate_variance(std::map<std::string, std::vector<float> > distances,map<std::string, std::vector<float> >& viriation);
void load_all_methylation_sites(const char* filename, std::map<std::string, std::vector<float> >& count_for_var);
std::vector<float> find_variance(std::vector<float> variance);
std::vector<float> find_variance(std::vector<float> variance);
void cmp_methylation(const char* filename, std::map<std::string, std::vector<float> > sites_mean, std::map<std::string, float>& aftermap);

bool sort_by_value(const std::pair<string,float>& a,const std::pair<string,float>& b){
	return a.second > b.second;
}


void sort_and_output_map(const char* filename, std::map<std::string, float>& aftermap, float upper) {
	std::ofstream outfile(filename);
	std::map<std::string, float>::iterator it;
	std::vector<std::pair<std::string, float> > sorted_map;
	for (it = aftermap.begin(); it != aftermap.end(); it++) {
		sorted_map.push_back(make_pair(it->first, it->second));
	}
	sort(sorted_map.begin(), sorted_map.end(), sort_by_value);
	for (std::vector<std::pair<std::string, float> >::iterator it2 = sorted_map.begin(); it2 != sorted_map.end(); it2++) {
		if (it2->second > upper) {
			if(it2->second < 1){
				outfile << it2->first << '\t' << it2->second << '\n';
			}
		}
		
	}
	std::cout << "Output complete" << endl;
	outfile.close();
}


void load_all_methylation_sites(const char* filename, std::map<std::string, std::vector<float> >& count_for_var) {
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
	while (*p != '\n' && (p - buffer) < size) { p++; }
	p++;
	while ((p - buffer) < size) {
		q = p;
		while (*p != '\t' && (p - buffer) < size) { p++; }
		string site_name = "";
		while (q < p) {
			site_name += *q;
			q++;
		}
		p++;
		q = p;
		if (*p != 'N') {
			while (*p != '\t' && *p != '\n' && (p - buffer) < size) { p++; }
			char *num = new char();
			memcpy(num, q, p - q);
			float beta = atof(num);
			delete[] num;
			count_for_var[site_name].push_back(beta);
		}
		while (*p != '\n' && (p - buffer) < size) { p++; }
		p++;
	}
	delete[] buffer;
	std::cout << "Map file load over." << endl;

}

std::vector<float> find_variance(std::vector<float> variance) {
	std::vector<float>::iterator it;
	float sum = 0;
	for (it = variance.begin(); it != variance.end(); it++) {
		sum += *it;
	}
	float mean = sum / variance.size();
	float stdev = 0;
	for (it = variance.begin(); it != variance.end(); it++) {
		stdev += (*it - mean) * (*it - mean);
	}
	stdev /= variance.size();
	std::vector<float> result;
	result.push_back(mean);
	result.push_back(sqrt(stdev));
	return result;
}


std::vector<float> generate_input_data(std::map<std::string, float> sites_mean) {
	std::vector<float> outputs;
	for (std::map<std::string, float>::iterator it = sites_mean.begin(); it != sites_mean.end(); it++) {
		if (it->second <= 0)outputs.push_back(1e-9);
		else if (it->second < 1)outputs.push_back(it->second);
	}
	return outputs;
}

void load_mean_vir_and_site(const char* filename, std::map<std::string, std::vector<float> >& mean_vir, std::map<std::string, int>& rgsites, int sites_size = 10000000) {
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
	int cow_number = 0;
	while ((p - buffer) < size) {
		std::vector<float> param;
		q = p;
		while (*p != '\t' && (p - buffer) < size) { p++; }
		std::string site_name = "";
		while (q<p) {
			site_name += *q;
			q++;
		}
		p++;
		q = p;
		while (*p != '\t' && (p - buffer) < size) { p++; }
		char *num = new char();
		memcpy(num, q, p - q);
		float mean = atof(num);
		delete[] num;
		param.push_back(mean);
		p++;
		q = p;
		while (*p != '\n' && (p - buffer) < size) { p++; }
		char *num2 = new char();
		memcpy(num2, q, p - q);
		float vir = atof(num2);
		delete[] num2;
		param.push_back(vir);
		p++;
		cow_number++;
		if(rgsites.size() != 0){
			if(rgsites.find(site_name) == rgsites.end()){
				mean_vir.insert(make_pair(site_name, param));
			}
		}
		else{
			mean_vir.insert(make_pair(site_name, param));
		}
		if (cow_number >= sites_size) {
			break;
		}
	}
	delete[] buffer;
	std::cout << "Map file load over." << endl;
}

void cmp_methylation(const char* filename, std::map<std::string, std::vector<float> > sites_mean, std::map<std::string, float>& aftermap) {
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
	std::string head_line = "";
	char *p = buffer;
	char* q = p;
	while (*p != '\n') { p++; }
	while (q < p) {
		head_line += *q;
		q++;
	}
	p++;
	q = p;
	/*std::string temp_head_line = "Composite Element REF	Beta_value	Chromosome	Start	End	Gene_Symbol	Gene_Type	Transcript_ID	Position_to_TSS	CGI_Coordinate	Feature_Type";
	if (head_line != temp_head_line) {
		std::cerr << "fatal error: The input file format is wrong, please check whether the correct input file is entered.\n";
		exit(1);
	}*/
	float beta;
	while ((p - buffer) < size) {
		q = p;
		while (*p != '\t' && (p - buffer) < size) { p++; }
		std::string site_name = "";
		while (q<p) {
			site_name += *q;
			q++;
		}
		p++;
		q = p;
		if (*p == '0' || *p == '1') {
			while (*p != '\t' && *p != '\n' && (p - buffer) < size) { p++; }
			char *num = new char();
			memcpy(num, q, p - q);
			beta = atof(num);
			delete[] num;
			std::map<std::string, std::vector<float> >::iterator it = sites_mean.find(site_name);
			if (it != sites_mean.end()) {
				if (sites_mean[site_name][0] < 0.5) {
					if ((std::abs(sites_mean[site_name][0] - beta) > 20 * sites_mean[site_name][1]))
						aftermap.insert(make_pair(site_name, abs(sites_mean[site_name][0] - beta) / (0.979689 - sites_mean[site_name][0])));
				}
				else {
					if ((std::abs(sites_mean[site_name][0] - beta) > 20 * sites_mean[site_name][1]))
						aftermap.insert(make_pair(site_name, abs(sites_mean[site_name][0] - beta) / (sites_mean[site_name][0] - 0.0219924)));
				}
			}
		}
		while (*p != '\n' && (p - buffer) < size) { p++; }
		p++;
	}
	delete[] buffer;
	std::cout << "Compair complete." << endl;
}


float calculate_purity(std::vector<float>& my_data, UsrParameters &UsrParameter){
	int data_size = my_data.size();
	float purity = 0;
	bmm_results bmm_result;
	if(data_size == 0){
	}
	else if(data_size == 1){
		purity = my_data[0];
	}
	else if(data_size >= 2 && data_size <= 3){
		kmeans_output ansmeans = my_kmeans_little(my_data);
		purity = ansmeans._centers[0];
	}
	else{ 
		if(data_size < 2 * UsrParameter.ClassNumber){
			UsrParameter.ClassNumber = data_size/2;
		}
		hyperparameter hyperparameters;
		init_bmm_hyperparameters(hyperparameters, UsrParameter.ClassNumber);
		parameter parameters;
		init_bmm_parameters(parameters, hyperparameters, my_data, UsrParameter.ClassNumber);
		bmm_result = bmm_function(my_data, UsrParameter.ClassNumber, parameters, hyperparameters, UsrParameter.RemoveCutoff, UsrParameter.MaxiterTime);
		if(UsrParameter.verbose){
			std::cout << "Bmm clusters number: " << bmm_result.cluster_num.size() << " ." << endl;
			std::cout << "Bmm clusters centers: ";
			std::cout << bmm_result.cluster_mean << endl;
		}
		purity = *max_element(bmm_result.cluster_mean.begin(),bmm_result.cluster_mean.end());
	}
	return purity;
}


std::vector<std::string> read_file_names(const char* filename){
	std::vector<std::string> file_names;
	std::ifstream fin(filename);
	std::string line;
	while (std::getline(fin, line)){
		file_names.push_back(line);
	}
	fin.close();
	return file_names;
}

void load_rgsites(const char* filename, std::map<std::string, int>& rgsites){
	filebuf *pbuf;
	ifstream filestr;
	long size;
	char *buffer;
	filestr.open(filename, ios::binary);
	//puts(filename);
	pbuf = filestr.rdbuf();
	size = pbuf->pubseekoff(0, ios::end, ios::in);
	pbuf->pubseekpos(0, ios::in);
	//char buffer[size];
	// 分配内存空间  
	
	buffer = new char[size];

	// 获取文件内容  
	pbuf->sgetn(buffer, size);

	filestr.close();
	char *p = buffer;
	char* q = p;
	while ((p - buffer)<size) {
		q = p;
		while (*p != '\n' && (p - buffer)<size) { p++; }
		std::string site_name = "";
		while (q<p) {
			site_name += *q;
			q++;
		}
		rgsites.insert(make_pair(site_name, 1));
		p++;
		q = p;
	}
	delete[] buffer;
}

#endif