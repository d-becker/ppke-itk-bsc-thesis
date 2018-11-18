// Computational loop 3 - Average density of each material over neighborhood of each cell
for (int j = 1; j < sizey-1; j++) {
    for (int i = 1; i < sizex-1; i++) {
	// o: outer
	double xo = ccc.x[i+sizex*j];
	double yo = ccc.y[i+sizex*j];

	// There are at most 9 neighbours in 2D case.
	double dsqr[9];

	// for all neighbours
	for (int nj = -1; nj <= 1; nj++) {
	    for (int ni = -1; ni <= 1; ni++) {

		dsqr[(nj+1)*3 + (ni+1)] = 0.0;

		// i: inner
		double xi = ccc.x[(i+ni)+sizex*(j+nj)];
		double yi = ccc.y[(i+ni)+sizex*(j+nj)];

		dsqr[(nj+1)*3 + (ni+1)] += (xo - xi) * (xo - xi);
		dsqr[(nj+1)*3 + (ni+1)] += (yo - yi) * (yo - yi);
	    }
	}

	int ix = ccc.imaterial[i+sizex*j];

	if (ix <= 0) {
	    // condition is 'ix >= 0', this is the equivalent of
	    // 'until ix < 0' from the paper
	    for (ix = -ix; ix >= 0; ix = ccc.nextfrac[ix]) {
		    int mat = ccc.matids[ix];
		    double rho_sum = 0.0;
		    int Nn = 0;

		    // for all neighbours
		    for (int nj = -1; nj <= 1; nj++) {
			for (int ni = -1; ni <= 1; ni++) {
			    int ci = i+ni, cj = j+nj;
			    int jx = ccc.imaterial[ci+sizex*cj];

			    if (jx <= 0) {
				// condition is 'jx >= 0', this is the equivalent of
				// 'until jx < 0' from the paper
				for (jx = -jx; jx >= 0; jx = ccc.nextfrac[jx]) {
					if (ccc.matids[jx] == mat) {
					    rho_sum += ccc.rho_compact_list[jx] /
						dsqr[(nj+1)*3 + (ni+1)];
					    Nn += 1;

					    // The loop has an extra condition: "and not found".
					    // This makes sense, if the material is found,
					    // there won't be any more of the same.
					    break;
					}
				    }
				}
				else {
				    // NOTE: In this case, the neighbour is a pure cell,
				    // its material index is in jx.
				    // NOTE: we index materials from zero,
				    // but zero can be a list index
				    int mat_neighbour = jx - 1;
				    if (mat == mat_neighbour) {
					rho_sum += ccc.rho_compact[ci+sizex*cj] /
					    dsqr[(nj+1)*3 + (ni+1)];
					Nn += 1;
				    }
				} // end if (jx <= 0)
			    } // end for (int ni)
			} // end for (int nj)

			ccc.rho_mat_ave_compact_list[ix] = rho_sum / Nn;
		    } // end for (ix = -ix)
		} // end if (ix <= 0)
		else {
		    // NOTE: In this case, the cell is a pure cell, its material index is in ix.
		    // NOTE: we index materials from zero, but zero can be a list index
		    int mat = ix - 1;

		    double rho_sum = 0.0;
		    int Nn = 0;

		    // for all neighbours
		    for (int nj = -1; nj <= 1; nj++) {
			if ((j + nj < 0) || (j + nj >= sizey))
			    continue;

			for (int ni = -1; ni <= 1; ni++) {
			    if ((i + ni < 0) || (i + ni >= sizex))
				continue;

			    int ci = i+ni, cj = j+nj;
			    int jx = ccc.imaterial[ci+sizex*cj];

			    if (jx <= 0) {
				// condition is 'jx >= 0', this is the equivalent of
				// 'until jx < 0' from the paper
				for (jx = -jx; jx >= 0; jx = ccc.nextfrac[jx]) {
					if (ccc.matids[jx] == mat) {
					    rho_sum += ccc.rho_compact_list[jx] /
						dsqr[(nj+1)*3 + (ni+1)];
					    Nn += 1;

					    // The loop has an extra condition: "and not found".
					    // This makes sense, if the material is found,
					    // there won't be any more of the same.
					    break;
					}
				    }
				}
				else {
				    // NOTE: In this case, the neighbour is a pure cell,
				    // its material index is in jx.
				    // NOTE: we index materials from zero,
				    // but zero can be a list index
				    int mat_neighbour = jx - 1;
				    if (mat == mat_neighbour) {
					rho_sum += ccc.rho_compact[ci+sizex*cj] /
					    dsqr[(nj+1)*3 + (ni+1)];
					Nn += 1;
				    }
				} // end if (jx <= 0)
			    } // end for (int ni)
			} // end for (int nj)

			ccc.rho_mat_ave_compact[i+sizex*j] = rho_sum / Nn;
		    } // end else
		}
	    }
	}
    }
}
