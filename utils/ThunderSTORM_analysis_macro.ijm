// Set path-parameters
points = "p50_sigma2"
points_per_gene = "50"
amplitude = "100"
bkg = "1000"

// Load image stack
run("Bio-Formats Importer", "open=C:/Users/jakob/Documents/Repositories/starfish_simulation_analysis/data/"+points+"/a"+amplitude+"_p"+points_per_gene+"_bkg"+bkg+"/primary/primary-fov_000-c0-r0-z0.tiff color_mode=Default group_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT dimensions axis_1_number_of_images=16 axis_1_axis_first_image=0 axis_1_axis_increment=1 contains=[] name=C:/Users/jakob/Documents/Repositories/starfish_simulation_analysis/data/"+points+"/a"+amplitude+"_p1_bkg"+bkg+"/primary/primary-fov_000-c0-r0-z0.tiff");
resetMinAndMax();

// Run ThunderSTORM analysis
run("Camera setup", "offset=0.0 isemgain=false photons2adu=3.6 pixelsize=1.0");
run("Run analysis", "filter=[Lowered Gaussian filter] sigma=1.6 detector=[Local maximum] connectivity=8-neighbourhood threshold=std(Wave.F1) estimator=[PSF: Integrated Gaussian] sigma=3.0 fitradius=3 method=[Maximum likelihood] full_image_fitting=false mfaenabled=false renderer=[No Renderer]");
run("Export results", "filepath=C:\\Users\\jakob\\Documents\\Repositories\\MERFISH_kNNDecoder\\Data\\"+points+"\\a"+amplitude+"_p"+points_per_gene+"_bkg"+bkg+"\\tsoutput.csv fileformat=[CSV (comma separated)] sigma=true intensity=true offset=true saveprotocol=false x=true y=true bkgstd=true id=true uncertainty=true frame=true");

