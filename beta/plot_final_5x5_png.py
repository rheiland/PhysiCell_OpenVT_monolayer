import matplotlib.pyplot as plt
import matplotlib.image as mpimg
# import numpy as np

# montage -geometry +0+0 -tile 5x5 out_cell_area_b0.09078_g0.93674/keep.png out_cell_area_b0.766_g0.93674/keep.png out_cell_area_b0.94701_g0.93674/keep.png out_cell_area_b0.97044_g0.93674/keep.png out_cell_area_b0.9866_g0.93674/keep.png out_cell_area_b0.09078_g0.9081/keep.png out_cell_area_b0.766_g0.9081/keep.png out_cell_area_b0.94701_g0.9081/keep.png out_cell_area_b0.97044_g0.9081/keep.png out_cell_area_b0.9866_g0.9081/keep.png out_cell_area_b0.09078_g0.77982/keep.png out_cell_area_b0.766_g0.77982/keep.png out_cell_area_b0.94701_g0.77982/keep.png out_cell_area_b0.97044_g0.77982/keep.png out_cell_area_b0.9866_g0.77982/keep.png out_cell_area_b0.09078_g0.68266/keep.png out_cell_area_b0.766_g0.68266/keep.png out_cell_area_b0.94701_g0.68266/keep.png out_cell_area_b0.97044_g0.68266/keep.png

# gamma_vals = [0.0, 0.68266, 0.77982, 0.9081, 0.93674]  # rows
# beta_vals = [0.09078, 0.766, 0.94701, 0.97044, 0.9866]  # columns
# gamma_vals = [0.0, 0.683, 0.816, 0.908, 0.937]
# beta_vals  = [0.0, 0.970, 0.987, 0.996, 0.998]

beta_vals  = [0.0, 0.978, 0.988, 0.994, 0.996]   # use 0.0 instead of 0.499 ??
gamma_vals = [0.0, 0.711, 0.826, 0.903, 0.911]
gamma_vals = [0.0, 0.711, 0.826]

gamma_vals.reverse()

# 1. Configuration
nrows, ncols = 5, 5
nrows, ncols = 3, 5
# Replace with your 25 image paths
# image_files = [f"image_{i}.png" for i in range(nrows * ncols)]

image_files = [f"out_cell_area_b{beta}_g{gamma}/keep.png" for gamma in gamma_vals for beta in beta_vals]
print(image_files)

# 2. Create subplots with shared axes
fig_w = 8
fig_h = 8
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_w, fig_h),
                        sharex=True, sharey=True)

# 3. Iterate through axes and images
for i in range(nrows):
# for i in [4,3,2,1,0]:
    gamma = gamma_vals[i]
    for j in range(ncols):
        beta = beta_vals[j]
        idx = i * ncols + j
        # print("idx=",idx)
        ax = axs[i, j]
        
        # Load and display image (assuming images exist)
        try:
            img = mpimg.imread(image_files[idx])
            ax.imshow(img)
        except FileNotFoundError:
            ax.text(0.5, 0.5, 'No Image', ha='center', va='center')
            
        # Hide ticks for inner plots, only show on edge
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        
        # Optional: Add row/column labels
        # if i == 0: ax.set_title(f"Col {j+1}")
        if j == 0: ax.set_ylabel(f"Row {i+1}")
        # if j == 0: ax.set_xlabel(r'$\beta$='+str(beta))
        if i==nrows-1 and j == 0: ax.set_xlabel(r'$\beta$='+str(beta))

# 4. Final adjustments
# plt.subplots_adjust(wspace=0.05, hspace=0.05) # Minimize gap
plt.savefig("final_5x5.pdf")   # , bbox_inches='tight', pad_inches=0)
plt.show()

