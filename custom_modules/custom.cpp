/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "custom.h"
#include "../BioFVM/BioFVM.h" 
using namespace BioFVM;

// declare cell definitions here 

std::vector<bool> nodes;

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = micropattern_update_cell_velocity;//NULL; // custom velocity function used, NULL for no motility
	cell_defaults.functions.update_phenotype = NULL; 
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.pre_update_intracellular = pre_update_intracellular; 
	cell_defaults.functions.post_update_intracellular = post_update_intracellular; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0 ); //for paraview visualization

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries();

	/*
	   This summarizes the setup. 
	*/
	
	display_cell_definitions( std::cout ); 


	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	int i = microenvironment.find_density_index("BMP4");
	int j = microenvironment.find_density_index("WNT");
	int k = microenvironment.find_density_index("NODAL");
	setup_densities(i,j,k);	// non-homogeneous initial conditions.
	return; 
}

void setup_tissue( void )
{
	// load cells from your CSV file
	load_cells_from_pugixml(); 	
}

void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	/*if (PhysiCell::PhysiCell_globals.current_time >= 100.0 
		&& pCell->phenotype.intracellular->get_parameter_value("$time_scale") == 0.0
	){
		pCell->phenotype.intracellular->set_parameter_value("$time_scale", 0.1);
	}*/
	static int WNT_index = microenvironment.find_density_index( "WNT" );
	static double WNT_threshold = parameters.doubles("WNT_threshold");
	static int BMP4_index = microenvironment.find_density_index( "BMP4" );
	static double BMP4_threshold = parameters.doubles("BMP4_threshold");
	static int NODAL_index = microenvironment.find_density_index( "NODAL" );
	static double NODAL_threshold = parameters.doubles("NODAL_threshold");
	
	static double WNT_saturation = parameters.doubles("WNT_saturation");
	static double BMP4_saturation = parameters.doubles("BMP4_saturation");
	static double NODAL_saturation = parameters.doubles("NODAL_saturation");
	
	static bool source_on = parameters.bools("source_when_on");
	static double WNT_source_amount = parameters.doubles("WNT_source_amount");
	static double BMP4_source_amount = parameters.doubles("BMP4_source_amount");
	static double NODAL_source_amount = parameters.doubles("NODAL_source_amount");

	if (WNT_index != -1)
	{
		double WNT_cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[WNT_index];
		pCell->phenotype.intracellular->set_boolean_variable_value("WNT", WNT_cell_concentration >= WNT_threshold);
		if (WNT_cell_concentration > WNT_saturation){pCell->phenotype.molecular.internalized_total_substrates[WNT_index] = WNT_saturation;}
		
		if (pCell->phenotype.intracellular->get_boolean_variable_value("WNT") && source_on){
			pCell->nearest_density_vector()[WNT_index] += WNT_source_amount;
		}
	} 
	if (BMP4_index != -1)
	{
		double BMP4_cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[BMP4_index];
		pCell->phenotype.intracellular->set_boolean_variable_value("BMP4", BMP4_cell_concentration >= BMP4_threshold);
		if (BMP4_cell_concentration > BMP4_saturation){pCell->phenotype.molecular.internalized_total_substrates[BMP4_index] = BMP4_saturation;}
		
		if (pCell->phenotype.intracellular->get_boolean_variable_value("BMP4") && source_on){
			pCell->nearest_density_vector()[BMP4_index] += BMP4_source_amount;
		}
	}
	
	if (NODAL_index != -1)
	{
		double NODAL_cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[NODAL_index];
		pCell->phenotype.intracellular->set_boolean_variable_value("NODAL", NODAL_cell_concentration >= NODAL_threshold);	
		if (NODAL_cell_concentration > NODAL_saturation){pCell->phenotype.molecular.internalized_total_substrates[NODAL_index] = NODAL_saturation;}
		
		if (pCell->phenotype.intracellular->get_boolean_variable_value("NODAL") && source_on){
			pCell->nearest_density_vector()[NODAL_index] += NODAL_source_amount;
		}
	}
	else{printf("NO INDICES \n");} // this would be an error, indices should be > 0.
}

void post_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	static double v_th = parameters.doubles("v_th");	
	static int NODAL_index = microenvironment.find_density_index( "NODAL" );
	double NODAL_cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[NODAL_index];	
	static double NODAL_threshold = parameters.doubles("NODAL_threshold");
	
	if (NODAL_cell_concentration >= v_th){
		pCell->phenotype.intracellular->set_parameter_value("$auto_inhib", true); // overides the logic to make NODAL autocatalytic		
	}
	else{
		pCell->phenotype.intracellular->set_parameter_value("$auto_inhib", false);
	}
	//printf("auto_inhib = %f \n", pCell->phenotype.intracellular->get_parameter_value("$auto_inhib"));
	
	color_node(pCell);
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
	if (parameters.bools("multi_node_colour")){ //see the Venn diagram in ./config for an overview of colour scheme. This func in not scaleable to larger GRNs
		if ( !pCell->phenotype.intracellular->get_boolean_variable_value( "BMP4" )  && !pCell->phenotype.intracellular->get_boolean_variable_value( "WNT" ) && !pCell->phenotype.intracellular->get_boolean_variable_value( "NODAL" )){
			output[0] = "rgb(125,125,125)"; // grey, no nodes are true
			output[2] = "rgb(100,100,100)";
		}
		
		if ( pCell->phenotype.intracellular->get_boolean_variable_value( "BMP4" )  && !pCell->phenotype.intracellular->get_boolean_variable_value( "WNT" ) && !pCell->phenotype.intracellular->get_boolean_variable_value( "NODAL" )){
			output[0] = "rgb(255,0,0)"; // red, only BMP4 = True
			output[2] = "rgb(125,0,0)";
		}
		if ( pCell->phenotype.intracellular->get_boolean_variable_value( "BMP4" )  && pCell->phenotype.intracellular->get_boolean_variable_value( "WNT" ) && !pCell->phenotype.intracellular->get_boolean_variable_value( "NODAL" )){
			output[0] = "rgb(255,255,0)"; // yellow, BMP4, WNT = True
			output[2] = "rgb(255,125,0)";
		}
		if ( pCell->phenotype.intracellular->get_boolean_variable_value( "BMP4" )  && !pCell->phenotype.intracellular->get_boolean_variable_value( "WNT" ) && pCell->phenotype.intracellular->get_boolean_variable_value( "NODAL" )){
			output[0] = "rgb(255,0,255)"; // purple, BMP4, NODAL = True
			output[2] = "rgb(125,0,125)";
		}
		
		if ( !pCell->phenotype.intracellular->get_boolean_variable_value( "BMP4" )  && pCell->phenotype.intracellular->get_boolean_variable_value( "WNT" ) && !pCell->phenotype.intracellular->get_boolean_variable_value( "NODAL" )){
			output[0] = "rgb(0,255,0)"; // green, only WNT = True
			output[2] = "rgb(0,125,0)";
		}
		if ( !pCell->phenotype.intracellular->get_boolean_variable_value( "BMP4" )  && pCell->phenotype.intracellular->get_boolean_variable_value( "WNT" ) && pCell->phenotype.intracellular->get_boolean_variable_value( "NODAL" )){
			output[0] = "rgb(0,255,255)"; // teal, only WNT, NODAL = True
			output[2] = "rgb(0,125,125)";
		}
		
		if ( !pCell->phenotype.intracellular->get_boolean_variable_value( "BMP4" )  && !pCell->phenotype.intracellular->get_boolean_variable_value( "WNT" ) && pCell->phenotype.intracellular->get_boolean_variable_value( "NODAL" )){
			output[0] = "rgb(0,0,255)"; // blue, only NODAL = True
			output[2] = "rgb(0,0,125)";
		}
		if ( pCell->phenotype.intracellular->get_boolean_variable_value( "BMP4" )  && pCell->phenotype.intracellular->get_boolean_variable_value( "WNT" ) && pCell->phenotype.intracellular->get_boolean_variable_value( "NODAL" )){
			output[0] = "rgb(90,30,10)"; // brown, all = True
			output[2] = "rgb(45,15,5)";
		}

	}
	else { //original colouring based on one node
		if ( !pCell->phenotype.intracellular->get_boolean_variable_value( parameters.strings("node_to_visualize") ) )
		{
			output[0] = "rgb(255,0,0)";
			output[2] = "rgb(125,0,0)";
		
		}
		else{
			output[0] = "rgb(0, 255,0)";
			output[2] = "rgb(0, 125,0)";
		}
	}
	
	return output;
}

void color_node(Cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	//pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}

void maintain_BMP4_density(int BMP4_index){ // this is deprecated in favour of maintain_signal_density
	static double r_loc = parameters.doubles("r_loc"); // see https://doi.org/10.1371/journal.pbio.3000498.s030
	static double b_ext = parameters.doubles("b_ext"); 
	#pragma omp parallel for
	for( int n=0; n < microenvironment.number_of_voxels() ; n++ ){
		auto current_voxel = microenvironment.voxels(n);
		double b = microenvironment.density_vector(n)[BMP4_index];
		std::vector<double> cent = {current_voxel.center[0], current_voxel.center[1], current_voxel.center[2]};
		if ((r_loc - norm(cent)) <= 0){
			microenvironment.density_vector(n)[BMP4_index] = b_ext; 		
		}
	}
}

void maintain_signal_density(int BMP4_index, int WNT_index, int NODAL_index){
	static double r_loc = parameters.doubles("r_loc"); // see https://doi.org/10.1371/journal.pbio.3000498.s030
	static double b_ext = parameters.doubles("b_ext"); 
	static double deg_rate = parameters.doubles("enhanced_deg_frac");
	#pragma omp parallel for
	for( int n=0; n < microenvironment.number_of_voxels() ; n++ ){
		auto current_voxel = microenvironment.voxels(n);
		double u = microenvironment.density_vector(n)[WNT_index];
		double v = microenvironment.density_vector(n)[NODAL_index];
		std::vector<double> cent = {current_voxel.center[0], current_voxel.center[1], current_voxel.center[2]};
		if ((r_loc - norm(cent)) < 0){
			microenvironment.density_vector(n)[BMP4_index] = b_ext; 
			microenvironment.density_vector(n)[WNT_index] =  deg_rate * u;	 // degredation on WNT and NODAL is enhanced outside of the colony
			microenvironment.density_vector(n)[NODAL_index] =  deg_rate * v;	
		}
	}
}

void setup_densities(int BMP4_index, int WNT_index, int NODAL_index){ // random intial concentration can be used, but behaves weirdly with cells.
	static double r_loc = parameters.doubles("r_loc");
	static double delta = parameters.doubles("delta"); // thickness of initial WNT, NODAL ring
	#pragma omp parallel for
	for( int n=0; n < microenvironment.number_of_voxels() ; n++ ){
		auto current_voxel = microenvironment.voxels(n);
		double b_rand = (double) rand()/RAND_MAX;//(rand()/RAND_MAX - 0.5) * 2; // random number [0,+1]
		std::vector<double> cent = {current_voxel.center[0], current_voxel.center[1], current_voxel.center[2]};
		if ((r_loc - norm(cent)) <= 0){
			microenvironment.density_vector(n)[BMP4_index] = 1.0;
		}
		if ((r_loc - norm(cent)) <= delta && (r_loc - norm(cent)) >= 10){// 10 units from edge, this could be added to XML to change without make.
			microenvironment.density_vector(n)[WNT_index] = 1.0; 	
			microenvironment.density_vector(n)[NODAL_index] = 1.0; 
			microenvironment.density_vector(n)[BMP4_index] = 0;//b_rand;
		}
		else{			
			double u_rand = (double) rand()/RAND_MAX;//(rand()/RAND_MAX - 0.5) * 2;
			double v_rand = (double) rand()/RAND_MAX;//(rand()/RAND_MAX - 0.5) * 2;
			microenvironment.density_vector(n)[BMP4_index] = 0;//b_rand;
			microenvironment.density_vector(n)[WNT_index] = 0;//u_rand; 	
			microenvironment.density_vector(n)[NODAL_index] = 0;//v_rand;
		}
	}
}

void micropattern_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt) //standard velocity function with a check to ensure cells remain inside micropatterns
{
	static double pattern_radius = parameters.doubles("r_loc"); //r_loc ~= pattern radius
	
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	pCell->state.neighbors.clear(); // new 1.8.0
	
	//First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		pCell->add_potentials(*neighbor);
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			pCell->add_potentials(*neighbor);
		}
	}
	
	//Main modifications here
	pCell->update_motility_vector(dt); 
	std::vector<double> vel_temp = pCell->velocity;
	pCell->velocity += phenotype.motility.motility_vector; 
	std::vector<double> vel_new = pCell->velocity;
	std::vector<double> cell_position = pCell->position;
	std::vector<double> dr = {vel_new[0]*dt, vel_new[1]*dt, vel_new[2]*dt};
	if((pattern_radius - norm(cell_position + dr)) < 0){
		pCell->velocity = -1 * vel_temp;
	}
	
	return; 
}