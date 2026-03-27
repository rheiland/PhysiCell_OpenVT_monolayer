
# to find the actual f_i,a_i values, knowing the chosen CDF percentiles: (from root dir)
python beta/chosen_CDF_pct.py 100 cdf_1000cells_linear_growth  bg00_cells1000_   f_i    # or “a_i” as last arg

 python param_sweep.py ../project    # _no_diffusion

- when done:
 python ../beta/plot_all_new_frames.py

- while running, from another terminal:
 python ../beta/plot_cell_scalars_4states.py -s beta_or_gamma --show_colorbar -o out_cell_area_b0.9866_g0.9081 -f -1


- but still results in a lousy res .pdf
 python ../beta/plot_final_5x5_png.py 
