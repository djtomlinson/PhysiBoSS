/*
 * custom.cpp
 *
 *  Created on: 15 feb. 2023
 *      Author: Arnau Montagud
 *  Description: 
*/

#include <sstream>

#include "../core/PhysiCell.h"
#include "../core/PhysiCell_constants.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

// setup functions to help us along 

void create_cell_types( void );
void setup_tissue( void ); 

// set up the BioFVM microenvironment 
void setup_microenvironment( void ); 
void change_dirichlet_nodes( void ); 

// helper function to read init files
std::vector<std::vector<double>>  read_cells_positions(std::string filename, char delimiter, bool header);

int voxel_index( int i, int j, int k );

void add_dirichlet_node( int voxel_index, std::vector<double>& value );

void set_substrate_dirichlet_activation( int substrate_index , bool new_value );

// custom pathology coloring function 

std::vector<std::string> regular_colors( Cell* pCell );

// custom functions can go here 

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt );
void custom_function( Cell* pCell, Phenotype& phenotype , double dt );

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt ); 

static const double EPSILON = std::numeric_limits<double>::epsilon();
void set_substrate_density ( void );
void add_ecm_interaction ( Cell* pC, int index_ecm, int index_voxel );
void custom_update_velocity ( Cell* pCell, Phenotype& phenotype, double dt );
double add_ecm_interaction_amadrid ( Cell* pCell, int index_ecm, int index_voxel );
void custom_update_velocity_amadrid ( Cell* pCell, Phenotype& phenotype, double dt );

double current_value( double min, double max, double percent );

std::vector<std::string> my_coloring_function( Cell* );
std::vector<std::string> ECM_coloring_function( Cell*);
std::vector<std::string> phase_coloring_function( Cell* );
std::vector<std::string> node_coloring_function( Cell* );
void SVG_plot_ecm( std::string filename , Microenvironment& M, double z_slice , double time, std::vector<std::string> (*cell_coloring_function)(Cell*), std::string sub );
