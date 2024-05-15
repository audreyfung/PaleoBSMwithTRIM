from .minerals_def import Minerals
# Create new mineral object
olivine = Minerals('Olivine', 4)

# paths and values
raw_data_path_array = ['/Users/szechingaudreyfung/paleodetector_with_trim/SRIM_raw_data/Olivine/O/', \
    '/Users/szechingaudreyfung/paleodetector_with_trim/SRIM_raw_data/Olivine/Si/', \
    '/Users/szechingaudreyfung/paleodetector_with_trim/SRIM_raw_data/Olivine/Fe/', \
    '/Users/szechingaudreyfung/paleodetector_with_trim/SRIM_raw_data/Olivine/Mg/']


derived_data_path_array = ['/Users/szechingaudreyfung/paleodetector_with_trim/SRIM_derived_data/Olivine/O/', \
    '/Users/szechingaudreyfung/paleodetector_with_trim/SRIM_derived_data/Olivine/Si/', \
    '/Users/szechingaudreyfung/paleodetector_with_trim/SRIM_derived_data/Olivine/Fe/', \
    '/Users/szechingaudreyfung/paleodetector_with_trim/SRIM_derived_data/Olivine/Mg/']

Er_avail_array = [0.01,0.03,0.05,0.07,0.1, 0.2, 0.3, 0.4, 0.5, \
                           5, 8, 10, 15, 18, 23, 25, 30, 33, 38, 40]

spin_array = [5/2, 1/2, 1/2, 5/2]


# set mineral attributes
olivine.set_composition(['O', 'Si', 'Fe', 'Mg'])
olivine.set_atomic_masses([15.999, 28.0855, 55.845, 24.305])
olivine.set_atomic_fractions([0.3863, 0.2674, 0.0547, 0.2916])
olivine.set_raw_data_paths(raw_data_path_array)
olivine.set_derived_data_paths(derived_data_path_array)
olivine.set_Er_avail(Er_avail_array)
olivine.set_nuclear_spins(spin_array)




print(olivine.name)