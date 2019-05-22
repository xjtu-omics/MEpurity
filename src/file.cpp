#include "bmm.h"
#include "help.h"
#include "kmeans.h"
#include "file.h"
using namespace std;

extern UsrParameters UsrParameter;




int main(int argc, char* argv[]) {
	Paramscan(argc, argv);
	std::map<std::string, std::vector<float> > sites_mean;
	std::map<std::string, int> rgsites;
	if(UsrParameter.rmSNPs != ""){
		load_rgsites(UsrParameter.rmSNPs.c_str(), rgsites);
	}
	load_mean_vir_and_site(UsrParameter.normalMapfile.c_str(), sites_mean, rgsites, UsrParameter.MapSitesNumber);
	std::ofstream outputt;
	outputt.open(UsrParameter.outputFilepath.c_str(), std::ios::app);
	if(UsrParameter.methylationfile != ""){
		std::map<string, float> dif_map;
		cmp_methylation(UsrParameter.methylationfile.c_str(), sites_mean, dif_map);
		std::vector<float> my_data = generate_input_data(dif_map);
		int data_size = my_data.size();
		float purity = calculate_purity(my_data, UsrParameter);
		outputt << UsrParameter.methylationfile << '\t' << purity << endl;
	}
	outputt.close();
	return 0;
}