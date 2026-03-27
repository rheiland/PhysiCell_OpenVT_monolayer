# from vis_base.py:

    def cell_counts_cb(self):
        # print("---- cell_counts_cb(): --> window for 2D population plots")
        # self.analysis_data_wait.value = 'compute n of N ...'

        if not self.get_cell_types_from_legend():
            if not self.get_cell_types_from_config():
                return

        xml_pattern = self.output_dir + "/" + "output*.xml"
        xml_files = glob.glob(xml_pattern)
        # print(xml_files)
        num_xml = len(xml_files)
        if num_xml == 0:
            print("last_plot_cb(): WARNING: no output*.xml files present")
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Information)
            msgBox.setText("Could not find any " + self.output_dir + "/output*.xml")
            msgBox.setStandardButtons(QMessageBox.Ok)
            msgBox.exec()
            return

        xml_files.sort()
        # print("sorted: ",xml_files)

        mcds = []
        for fname in xml_files:
            basename = os.path.basename(fname)
            # print("basename= ",basename)
            # mcds = pyMCDS(basename, self.output_dir, microenv=False, graph=False, verbose=False)
            mcds.append(pyMCDS(basename, self.output_dir, microenv=False, graph=False, verbose=False))

        if self.discrete_scalar not in mcds[0].data['discrete_cells']['data'].keys():
            print(f"\ncell_counts_cb(): {self.discrete_scalar} is not saved in the output. See the Full list above. Exiting.")
            return

        tval = np.linspace(0, mcds[-1].get_time(), len(xml_files))
        # print("  max tval=",tval)

        # self.yval4 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 4) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

        #--------
        if self.discrete_scalar == 'cell_type':   # number not known until run time
            # if not self.population_plot[self.discrete_scalar]:
            self.population_plot[self.discrete_scalar] = PopulationPlotWindow() # don't test if already exists!
            self.population_plot[self.discrete_scalar].ax0.cla()

            # ctype_plot = []
            lw = 2
            # for itype, ctname in enumerate(self.celltypes_list):
            # print("  self.celltype_name=",self.celltype_name)
            for itype in range(len(self.celltype_name)):
                ctname = self.celltype_name[itype]
                try:
                    ctcolor = self.celltype_color[itype]
                    # print("  cell_counts_cb(): ctcolor (1)=",ctcolor)
                    if ctcolor == "yellow":   # can't see yellow on white
                        ctcolor = "gold"
                        # print("  now cell_counts_cb(): ctcolor (1)=",ctcolor)
                except:
                    ctcolor = 'C' + str(itype)   # use random colors from matplotlib; TODO: avoid yellow, etc
                    # print("  cell_counts_cb(): ctcolor (2)=",ctcolor)
                if 'rgb' in ctcolor:
                    rgb = ctcolor.replace('rgb','')
                    rgb = rgb.replace('(','')
                    rgb = rgb.replace(')','')
                    rgb = rgb.split(',')
                    # print("--- rgb after split=",rgb)
                    ctcolor = [float(rgb[0])/255., float(rgb[1])/255., float(rgb[2])/255.]
                    # print("--- converted rgb=",ctcolor)
                if self.celltype_filter:
                    yval = np.array( [(np.count_nonzero((np.isin(mcds[idx].data['discrete_cells']['data']['cell_type'], self.celltype_filter)) & (mcds[idx].data['discrete_cells']['data']['cell_type'] == itype) & (mcds[idx].data['discrete_cells']['data']['cycle_model'] < 100.) == True)) for idx in range(len(mcds))] )
                else:
                    yval = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['data']['cell_type'] == itype) & (mcds[idx].data['discrete_cells']['data']['cycle_model'] < 100.) == True)) for idx in range(len(mcds))] )
                # yval = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['data']['cell_type'] == itype) == True)) for idx in range(len(mcds))] )
                # print("  yval=",yval)
                if yval.sum() > 0: # only plot if there are cells of this type
                    self.population_plot[self.discrete_scalar].ax0.plot(tval, yval, label=ctname, linewidth=lw, color=ctcolor)


            self.population_plot[self.discrete_scalar].ax0.set_xlabel('time (mins)')
            self.population_plot[self.discrete_scalar].ax0.set_ylabel('# of cells')
            self.population_plot[self.discrete_scalar].ax0.set_title("cell_type", fontsize=10)
            self.population_plot[self.discrete_scalar].ax0.legend(loc='center left', prop={'size': 8})
            self.population_plot[self.discrete_scalar].canvas.update()
            self.population_plot[self.discrete_scalar].canvas.draw()
            self.population_plot[self.discrete_scalar].show()

        #--------
        elif self.discrete_scalar == '"number_of_nuclei"':   # is it used yet?
            pass
        #--------
        else:  # number is fixed for these (cycle_model, current_phase, is_motile, current_death_model, dead)

            # [‘cell_type’, ‘cycle_model’, ‘current_phase’,‘is_motile’,‘current_death_model’,‘dead’,‘number_of_nuclei’,‘polarity’]
            # self.discrete_scalar_len = {"cell_type":0, "cycle_model":6, "current_phase":4, "is_motile":2,"current_death_model":2, "dead":2, "number_of_nuclei":0 }

            self.population_plot[self.discrete_scalar] = PopulationPlotWindow()
            self.population_plot[self.discrete_scalar].ax0.cla()

            # print("---- generate plot for ",self.discrete_scalar)
            # ctype_plot = []
            lw = 2
            # for itype, ctname in enumerate(self.celltypes_list):
            # print("  self.celltype_name=",self.celltype_name)
            # for itype in range(self.discrete_scalar_len[self.discrete_scalar]):
            for itype in self.discrete_scalar_vals[self.discrete_scalar]:
                # print("  cell_counts_cb(): itype= ",itype)
                ctcolor = 'C' + str(itype)   # use random colors from matplotlib
                # print("  ctcolor=",ctcolor)
                # yval = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['data']['cell_type'] == itype) & (mcds[idx].data['discrete_cells']['data']['cycle_model'] < 100.) == True)) for idx in range(len(mcds))] )

                # yval = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['data'][self.discrete_scalar] == itype) ) for idx in range(len(mcds))) ] )

                # yval = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['data'][self.discrete_scalar] == itype) & (mcds[idx].data['discrete_cells']['data']['cycle_model'] < 100.) == True)) for idx in range(len(mcds))] )
                # yval = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['data'][self.discrete_scalar] == itype) ) for idx in range(len(mcds)))] )
                # yval = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['data'][self.discrete_scalar] == itype) & True) for idx in range(len(mcds)))] )

                # TODO: fix this hackiness. Do we want to avoid counting dead cells??
                if self.discrete_scalar == 'current_death_model': # Hack: because current_death_model is not working in PhysiCell, using cycle_model instead  
                    if self.celltype_filter: # Cell type filter applied here
                        yval = np.array( [(np.count_nonzero((np.isin(mcds[idx].data['discrete_cells']['data']['cell_type'], self.celltype_filter)) & (mcds[idx].data['discrete_cells']['data']['cycle_model'] == itype) & (mcds[idx].data['discrete_cells']['data']['cycle_model'] < 999.) == True)) for idx in range(len(mcds))] )
                    else:
                        yval = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['data']['cycle_model'] == itype) & (mcds[idx].data['discrete_cells']['data']['cycle_model'] < 999.) == True)) for idx in range(len(mcds))] )
                else:
                    if self.celltype_filter: # Cell type filter applied here
                        yval = np.array( [(np.count_nonzero((np.isin(mcds[idx].data['discrete_cells']['data']['cell_type'], self.celltype_filter)) & (mcds[idx].data['discrete_cells']['data'][self.discrete_scalar] == itype) & (mcds[idx].data['discrete_cells']['data']['cycle_model'] < 999.) == True)) for idx in range(len(mcds))] )
                    else:
                        yval = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['data'][self.discrete_scalar] == itype) & (mcds[idx].data['discrete_cells']['data']['cycle_model'] < 999.) == True)) for idx in range(len(mcds))] )
                # print("  yval=",yval)

                # if (self.discrete_scalar == 'cycle_model'): mylabel = 
                # else:
                # Check if exist any cells in the entire simulation with self.discrete_scalar occuring
                mylabel = str(itype)
                bool_list = ['is_motile', 'dead']
                if( yval.sum() > 0 or self.discrete_scalar in bool_list): # only plot if there are cells with this scalar or boolean
                    if (self.discrete_scalar == 'cycle_model' or self.discrete_scalar == 'current_death_model'): mylabel = self.cycle_models[itype]
                    elif (self.discrete_scalar == 'current_phase'): mylabel = self.cycle_phases[itype]
                    elif (self.discrete_scalar in bool_list ): mylabel = str(bool(itype))
                    # Plot only if there are cells with this scalar
                    self.population_plot[self.discrete_scalar].ax0.plot(tval, yval, label=mylabel, linewidth=lw, color=ctcolor)
                # self.population_plot[self.discrete_scalar].ax0.plot(tval, yval, linewidth=lw, color=ctcolor)
                # print(self.discrete_scalar, itype, mylabel, yval.sum() )

            
            self.population_plot[self.discrete_scalar].ax0.set_xlabel('time (mins)')
            self.population_plot[self.discrete_scalar].ax0.set_ylabel('# of cells')
            if self.celltype_filter:
                self.population_plot[self.discrete_scalar].ax0.set_title(self.discrete_scalar + " (filtered by cell type)", fontsize=10)
            else:
                self.population_plot[self.discrete_scalar].ax0.set_title(self.discrete_scalar, fontsize=10)
            self.population_plot[self.discrete_scalar].ax0.legend(loc='center left', prop={'size': 8})
            self.population_plot[self.discrete_scalar].canvas.update()
            self.population_plot[self.discrete_scalar].canvas.draw()
            self.population_plot[self.discrete_scalar].show()
