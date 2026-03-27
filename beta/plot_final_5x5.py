import matplotlib.pyplot as plt
import matplotlib.image as mpimg
# import numpy as np

# montage -geometry +0+0 -tile 5x5 out_cell_area_b0.09078_g0.93674/keep.png out_cell_area_b0.766_g0.93674/keep.png out_cell_area_b0.94701_g0.93674/keep.png out_cell_area_b0.97044_g0.93674/keep.png out_cell_area_b0.9866_g0.93674/keep.png out_cell_area_b0.09078_g0.9081/keep.png out_cell_area_b0.766_g0.9081/keep.png out_cell_area_b0.94701_g0.9081/keep.png out_cell_area_b0.97044_g0.9081/keep.png out_cell_area_b0.9866_g0.9081/keep.png out_cell_area_b0.09078_g0.77982/keep.png out_cell_area_b0.766_g0.77982/keep.png out_cell_area_b0.94701_g0.77982/keep.png out_cell_area_b0.97044_g0.77982/keep.png out_cell_area_b0.9866_g0.77982/keep.png out_cell_area_b0.09078_g0.68266/keep.png out_cell_area_b0.766_g0.68266/keep.png out_cell_area_b0.94701_g0.68266/keep.png out_cell_area_b0.97044_g0.68266/keep.png

# gamma_vals = [0.0, 0.68266, 0.77982, 0.9081, 0.93674]  # rows
# beta_vals = [0.09078, 0.766, 0.94701, 0.97044, 0.9866]  # columns

beta_vals  = [0.0, 0.978, 0.988, 0.994, 0.996]   # use 0.0 instead of 0.499 ??
gamma_vals = [0.0, 0.711, 0.826, 0.903, 0.911]
gamma_vals = [0.0, 0.711, 0.826]

# 1. Create dummy image data (replace with your PNG loading logic)
# In reality, you'd load images using: image = plt.imread('path.png')
rows, cols = 5, 5
image_files = [f"image_{i}.png" for i in range(nrows * ncols)]


len_fig=8
fig, axs = plt.subplots(rows, cols, figsize=(len_fig, len_fig), sharex=True, sharey=True)

# 2. Loop through axes and display images
for i in range(rows):
    for j in range(cols):
        # Generate dummy data for demonstration
        # dummy_img = np.random.rand(100, 100) 
        
        # Plot the image
        axs[i, j].imshow(dummy_img, cmap='gray')
        
        # Remove tick labels for inner plots, keep on outer edges
        if i == rows - 1:
            axs[i, j].set_xlabel(f'Col {j+1}')
        if j == 0:
            axs[i, j].set_ylabel(f'Row {i+1}')
        
        # Hide tick marks for cleaner layout
        axs[i, j].set_xticks([])
        axs[i, j].set_yticks([])

# 3. Adjust spacing between subplots
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.suptitle('5x5 Image Grid', fontsize=16)
plt.savefig("final_5x5.pdf")   # , bbox_inches='tight', pad_inches=0)
plt.show()

