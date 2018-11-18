const std::size_t COLS = ...;
const std::size_t ROWS = ...;

// Specify which materials are present in each cell.
const std::vector<std::vector<std::size_t>> materials = ...;

Data<2> data({COLS, ROWS}, materials);

CellData<2> x = data.new_cell_data();
CellData<2> y = data.new_cell_data();
CellMataData<2> Vf = data.new_cell_mat_data();
CellMataData<2> rho = data.new_cell_mat_data();
CellMataData<2> rho_mat_ave = data.new_cell_mat_data();

// Fill the datasets with data.

void kernel(NeighProxy<CellData<2>> x,
	    NeighProxy<CellData<2>> y,
	    NeighProxy<CellMatData<2>> Vf,
	    NeighProxy<CellMatData<2>> rho,
	    double &rho_out) {
    double xo = x[{0,0}];
    double yo = y[{0,0}];

    double dsqr[9];

    for (int nj = -1; nj <= 1; nj++) {
	for (int ni = -1; ni <= 1; ni++) {

	    dsqr[(nj+1)*3 + (ni+1)] = 0.0;

	    double xi = x[{ni,nj}];
	    double yi = y[{ni,nj}];

	    dsqr[(nj+1)*3 + (ni+1)] += (xo - xi) * (xo - xi);
	    dsqr[(nj+1)*3 + (ni+1)] += (yo - yi) * (yo - yi);
	}
    }

    if (Vf[{0,0}] > 0.0) {
	double rho_sum = 0.0;
	int Nn = 0;

	for (int nj = -1; nj <= 1; nj++) {
	    for (int ni = -1; ni <= 1; ni++) {

		if (Vf[{ni,nj}] > 0.0) {
		    rho_sum += rho[{ni,nj}] / dsqr[(nj+1)*3 + (ni+1)];
		    Nn += 1;
		}
	    }
	}
	rho_out = rho_sum / Nn;
    } else {
	rho_out = 0.0;
    }
}

Stencil<2> s9pt({{1,1},  {1,0},  {1,-1},
		 {0,1},  {0,0},  {0,-1},
		 {-1,1}, {-1,0}, {-1,-1}});
IndexGenerator<2> index_generator({1, 1}, {sizex-1, sizey-1});
Computation<2> computation(data, index_generator);

computation.compute(kernel,
		    NEIGH<CellData<2>>(x, s9pt),
		    NEIGH<CellData<2>>(y, s9pt),
		    NEIGH<CellMatData<2>>(Vf, s9pt),
		    NEIGH<CellMatData<2>>(rho, s9pt),
		    OUT<CellMatData<2>>(rho_mat_ave));
