# Write code to make Euclid cutouts based on a given coordinate and saving them to a specified folder.
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.esa.euclid import EuclidClass

# Get access to the Euclid archive
Euclid = EuclidClass(environment='PDR')
Euclid.login(credentials_file='/data2/visserscott/Documents/Euclid/Euclid_tools/cred.txt') 

# Start with the functions required to get the files for the cutouts
def checkInput(instrument_name, nisp_filters):
    valid_instruments = ["VIS", "NISP"]
    valid_filters = ["NIR_H", "NIR_J", "NIR_Y", ""]
    if instrument_name not in valid_instruments:
        raise ValueError(f"Invalid instrument name. Choose from: {valid_instruments}") 
    for filter_name in nisp_filters:
        if filter_name not in valid_filters:
            raise ValueError(f"Invalid filter name. Choose from: {valid_filters}")     
    return
            
def chooseBestFile(df_res_dup):
    #choosing best file for sources with multiple (with pandas)
    min_dist_per_plotcode = df_res_dup.groupby(['right_ascension', 'filter_name'], as_index=False, sort=False)['dist'].idxmin()
    df_final = df_res_dup.loc[min_dist_per_plotcode["dist"]]
    return df_final


def getFiles(product_type, sources, instrument_name, nisp_filters=[], segmentation_map_id = False):
    #checking input value correctness
    checkInput(instrument_name, nisp_filters)
    
    if product_type == "mosaic":
        if segmentation_map_id == True:
            df = getFilesMosaic(sources, instrument_name, nisp_filters)
        else:
            df = getFilesMosaicBackup(sources, instrument_name, nisp_filters)        
    elif product_type == "calibrated_frame":
        df = getFilesCalib(sources, instrument_name, nisp_filters)
    else:
        raise ValueError(f"Invalid instrument name. Choose from: mosaic, calibrated_frame") 
        
    return df
    
def getFilesCalib(source_query, instrument_name, nisp_filters):
    
    query= "SELECT sources.*, frames.file_name, frames.file_path, frames.datalabs_path, frames.instrument_name, frames.filter_name, frames.observation_id, frames.ra, frames.dec, DISTANCE(frames.ra, frames.dec, sources.right_ascension, sources.declination) AS dist \
    FROM "+source_query+" as sources \
    JOIN sedm.calibrated_frame AS frames \
    ON (frames.instrument_name='"+instrument_name+"') AND (frames.fov IS NOT NULL AND CONTAINS(CIRCLE('ICRS', sources.right_ascension, sources.declination, 0.00833333), frames.fov)=1)"

    #modify query to include only the provided nisp filters
    if instrument_name == 'NISP' and len(nisp_filters) > 0 and len(nisp_filters) < 3:
        query = addFilters("frames", query, instrument_name, nisp_filters)
        
    #make query to get file info
    Euclid.login(credentials_file = '/data2/visserscott/Documents/Euclid/Euclid_tools/cred.txt')
    job = Euclid.launch_job_async(query, verbose=False)
    df_res = job.get_results().to_pandas()
    print("Found", len(df_res), "query results")
    
    return df_res
    
def getFilesMosaic(source_query, instrument_name, nisp_filters):
    
    query= "SELECT sources.*, mosaics.file_name, mosaics.file_path, mosaics.datalabs_path, mosaics.mosaic_product_oid, mosaics.tile_index, mosaics.instrument_name, mosaics.filter_name, mosaics.ra, mosaics.dec, DISTANCE(mosaics.ra, mosaics.dec, sources.right_ascension, sources.declination) AS dist \
        FROM "+source_query+" as sources \
        JOIN sedm.mosaic_product AS mosaics \
        ON (mosaics.instrument_name='"+instrument_name+"') AND (CAST(mosaics.tile_index AS CHAR(25)) = SUBSTRING(CAST(segmentation_map_id AS CHAR(25)), 1, 9))"

    #modify query to include only the provided nisp filters
    if instrument_name == 'NISP' and len(nisp_filters) > 0 and len(nisp_filters) < 3:
        query = addFilters("mosaics", query, instrument_name, nisp_filters)
        
    #make query to get file info
    Euclid.login(credentials_file = '/data2/visserscott/Documents/Euclid/Euclid_tools/cred.txt')
    job = Euclid.launch_job_async(query, verbose=False)
    
    if job == None:
        raise AttributeError("Query for files failed. Please make sure you are logged in, your source query has brackets around it and contains columns named 'right_ascension' and 'declination' ")
    
    df_res = job.get_results().to_pandas()
    print("Found", len(df_res), "query results")
    
    return df_res


def getFilesMosaicBackup(source_query, instrument_name, nisp_filters):
    
    query= "SELECT sources.*, mosaics.file_name, mosaics.file_path, mosaics.datalabs_path, mosaics.mosaic_product_oid, mosaics.tile_index, mosaics.instrument_name, mosaics.filter_name, mosaics.ra, mosaics.dec, DISTANCE(mosaics.ra, mosaics.dec, sources.right_ascension, sources.declination) AS dist \
        FROM "+source_query+" as sources \
        JOIN sedm.mosaic_product AS mosaics \
        ON (mosaics.instrument_name='"+instrument_name+"') AND (mosaics.fov IS NOT NULL AND CONTAINS(CIRCLE('ICRS', sources.right_ascension, sources.declination, 0.00833333), mosaics.fov)=1)"

    #modify query to include only the provided nisp filters
    if instrument_name == 'NISP' and len(nisp_filters) > 0 and len(nisp_filters) < 3:
        query = addFilters("mosaics", query, instrument_name, nisp_filters)
    
    #make query to get file info
    Euclid.login(credentials_file = '/data2/visserscott/Documents/Euclid/Euclid_tools/cred.txt')
    job = Euclid.launch_job_async(query, verbose=False)
    
    if job == None:
        raise AttributeError("Query for files failed. Please make sure you are logged in, your source query has brackets around it and contains columns named 'right_ascension' and 'declination' ")
    
    df_res_dup = job.get_results().to_pandas()
    df_res = chooseBestFile(df_res_dup)
    print("Found", len(df_res), "query results")
    
    return df_res


def addFilters(product_type, query, instrument_name, nisp_filters=[]):
    filter_part = ""
    query_parts = query.split("ON ")
    for one_filter in nisp_filters:
        if one_filter != nisp_filters[-1]:
            filter_part += "("+ product_type + ".filter_name='" + one_filter + "') OR "
        else:
            filter_part += "(" + product_type + ".filter_name='" + one_filter + "')"
    query = query_parts[0] + "ON (" + filter_part + ") AND" + query_parts[1]
    return query


# Now the main function to make the cutouts based on a given coordinate and saving them to a specified folder.
def make_euclid_cutout_with_coord(coord, instrument, cutout_radius, save_folder_path, nisp_filters=None):
    """
    Makes a cutout from the Euclid archive at the given coordinate.

    Parameters:
    -----------
    coord : astropy.coordinates.SkyCoord
        The sky coordinate where the cutout is centered.
    instrument : str
        The instrument to use for the cutout. Either "VIS" or "NISP" or "VIS + NISP".
    cutout_radius : astropy.units.Quantity
        The radius of the cutout (e.g., 2*u.arcmin).
    save_folder_path : str
        The folder path where the cutout will be saved.
    nisp_filters : list of str, optional
        List of NISP filters to use if instrument is "NISP". Default is None.

    Returns:
    --------
    list of str
        List of file paths to the saved cutouts.
    """
    # Check if coord is a SkyCoord object
    if not isinstance(coord, SkyCoord):
        raise ValueError("coord must be an astropy.coordinates.SkyCoord object")
    
    # Check if instrument is valid
    if instrument not in ["VIS", "NISP", "VIS + NISP"]:
        raise ValueError("instrument must be either 'VIS' or 'NISP' or 'VIS + NISP'")
    
    # Check if nisp_filters is provided when instrument is NISP
    if instrument in ["NISP", "VIS + NISP"] and (nisp_filters is None or not isinstance(nisp_filters, list)):
        raise ValueError("nisp_filters must be provided as a list when instrument is 'NISP' or 'VIS + NISP'")

    query = f"""(SELECT TOP 1 euclid.right_ascension, euclid.declination, euclid.object_id, euclid.segmentation_area, euclid.flux_vis_1fwhm_aper, euclid.segmentation_map_id,
            DISTANCE(POINT(euclid.right_ascension, euclid.declination), POINT({coord.ra.value}, {coord.dec.value}))*3600 as dist_arcsec
            FROM catalogue.mer_catalogue as euclid
            WHERE DISTANCE(POINT(euclid.right_ascension, euclid.declination), POINT({coord.ra.value}, {coord.dec.value})) < 10./3600)"""
    cutout_paths = []
    
    if instrument == "VIS":
        file_info = getFiles("mosaic", query, instrument, segmentation_map_id=True)
        filepath = file_info['file_path'].values[0]+"/"+file_info['file_name'].values[0]
        saved_cutout_filepath = Euclid.get_cutout(file_path = filepath, instrument = "VIS", id = file_info['segmentation_map_id'].values[0], 
                                                  coordinate = coord, radius = cutout_radius, 
                                                  output_file = f'{save_folder_path}Euclid_{instrument}_RA{coord.ra.deg}_DEC{coord.dec.deg}_cutout.fits')
        cutout_paths.append(saved_cutout_filepath)
        print(f"Finished making cutouts for coordinate RA: {coord.ra.deg}, DEC: {coord.dec.deg} for instrument: {instrument} with radius {cutout_radius}.")
    elif instrument == "NISP":
        file_info = getFiles("mosaic", query, instrument, nisp_filters=nisp_filters, segmentation_map_id=True)
        filepaths = [file_info['file_path'].values[i]+"/"+file_info['file_name'].values[i] for i in range(len(file_info))]
        for n in range(len(filepaths)):
            saved_cutout_filepath = Euclid.get_cutout(file_path = filepaths[n], instrument = "NISP", id = file_info['segmentation_map_id'].values[n], 
                                                      coordinate = coord, radius = cutout_radius, 
                                                      output_file = f'{save_folder_path}Euclid_{instrument}_{file_info["filter_name"].values[n]}_RA{coord.ra.deg}_DEC{coord.dec.deg}_cutout.fits')
            cutout_paths.append(saved_cutout_filepath)
            print(f"Finished making cutouts for coordinate RA: {coord.ra.deg}, DEC: {coord.dec.deg} for instrument: {instrument} for the filters {nisp_filters} with radius {cutout_radius}.")
    elif instrument == "VIS + NISP":
        # Create the VIS cutout
        file_info_vis = getFiles("mosaic", query, "VIS", segmentation_map_id=True)
        filepath_vis = file_info_vis['file_path'].values[0]+"/"+file_info_vis['file_name'].values[0]
        saved_cutout_filepath_vis = Euclid.get_cutout(file_path = filepath_vis, instrument = "VIS", id = file_info_vis['segmentation_map_id'].values[0], 
                                                  coordinate = coord, radius = cutout_radius, 
                                                  output_file = f'{save_folder_path}Euclid_VIS_RA{coord.ra.deg}_DEC{coord.dec.deg}_cutout.fits')
        cutout_paths.append(saved_cutout_filepath_vis)

        # Create the NISP cutouts
        file_info_nir = getFiles("mosaic", query, "NISP", nisp_filters=nisp_filters, segmentation_map_id=True)
        filepaths_nir = [file_info_nir['file_path'].values[i]+"/"+file_info_nir['file_name'].values[i] for i in range(len(file_info_nir))]
        for n in range(len(filepaths_nir)):
            saved_cutout_filepath_nir = Euclid.get_cutout(file_path = filepaths_nir[n], instrument = "NISP", id = file_info_nir['segmentation_map_id'].values[n], 
                                                      coordinate = coord, radius = cutout_radius, 
                                                      output_file = f'{save_folder_path}Euclid_NISP_{file_info_nir["filter_name"].values[n]}_RA{coord.ra.deg}_DEC{coord.dec.deg}_cutout.fits')
            cutout_paths.append(saved_cutout_filepath_nir)
        print(f"Finished making cutouts for coordinate RA: {coord.ra.deg}, DEC: {coord.dec.deg} for instrument: {instrument} for the NISP filters {nisp_filters} with radius {cutout_radius}.")

    return cutout_paths


##### For testing the function above #####
"""
# Define a test coordinate (RA, DEC in degrees)
test_ra = 266.41472093639214 # degree
test_dec = 66.25481639101852 # degree
test_coord = SkyCoord(ra=test_ra*u.degree, dec=test_dec*u.degree, frame='icrs')

# Example usage of the function use above for the test coordinate defined earlier
cutout_path_vis = make_euclid_cutout_with_coord(coord=test_coord, instrument="VIS", cutout_radius=2*u.arcmin, 
                                                save_folder_path='/data2/visserscott/Documents/Euclid/Euclid_tools/test_cutouts/')
cutout_path_nir_1filter = make_euclid_cutout_with_coord(coord=test_coord, instrument="NISP", cutout_radius=2*u.arcmin, 
                                                        save_folder_path='/data2/visserscott/Documents/Euclid/Euclid_tools/test_cutouts/', nisp_filters=['NIR_H'])
cutout_path_nir = make_euclid_cutout_with_coord(coord=test_coord, instrument="NISP", cutout_radius=2*u.arcmin, 
                                                save_folder_path='/data2/visserscott/Documents/Euclid/Euclid_tools/test_cutouts/', nisp_filters=['NIR_H', 'NIR_J', 'NIR_Y'])
cutout_path_both = make_euclid_cutout_with_coord(coord=test_coord, instrument="VIS + NISP", cutout_radius=2*u.arcmin, 
                                                 save_folder_path='/data2/visserscott/Documents/Euclid/Euclid_tools/test_cutouts/', nisp_filters=['NIR_H', 'NIR_J', 'NIR_Y'])
"""