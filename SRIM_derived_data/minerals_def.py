class Minerals:
    def __init__(self, mineral_name, number_of_elements):
        self.name = mineral_name
        self.number_of_elements = number_of_elements

    def set_number_of_elements(self, number):
        self.number_of_elements = number
    def set_composition(self, composition_array):
        """
        input should be a list of strings, say for Galena, the input should be:
        ['Pb', 'S']
        """
        if self.number_of_elements != len(composition_array):
            print('the length of you array does not match with the number of elements')
            return
        else:
            self.composition = composition_array
    
    def set_atomic_masses(self, atomic_mass_array):
        if self.number_of_elements != len(atomic_mass_array):
            print('the number of fractions you enter does not match with the number of elements')
            return 
        else:
            self.atomic_masses = atomic_mass_array

    def set_atomic_fractions(self, fraction_array):
        if self.number_of_elements != len(fraction_array):
            print('the number of fractions you enter does not match with the number of elements')
            return
        else:   
            self.atomic_fractions = fraction_array
    def set_nuclear_spins(self, spin_array):
        if self.number_of_elements != len(spin_array):
            print('the number of fractions you enter does not match with the number of elements')
            return
        else:   
            self.nuclear_spin = spin_array

    


    def set_raw_data_paths(self, raw_data_path_array):
        if self.number_of_elements != len(raw_data_path_array):
            print('the length of your array does not match with the number of elements')
            return
        else:
            self.raw_data_path = raw_data_path_array

    def set_derived_data_paths(self, derived_data_path_array):
        if self.number_of_elements != len(derived_data_path_array):
            print('the length of your array does not match with the number of elements')
            return 
        else:
            self.derived_data_path = derived_data_path_array

    def set_Er_avail(self, Er_avail_array):
        """
        input is a list of float of recoil energies where data files are available,
        here assumed all elements have the same avilable energies, so only accept 
        one array.
        """
        
        self.Er_avail = Er_avail_array

    
    def save_bezier_xt(self):
        """
        read in raw data and save data as standardized numpy array to be 
        loaded by load_bezier_xt
        """
        print('now in save_bezier_xt')

    def save_PCA_xt(self):
        """
        read in raw data and save data as standardized numpy array to be 
        loaded by load_PCA_xt
        """

    def load_bezier_xt(self):
        """
        go to the corresponding data paths to load the saved track length
        distributions calculated using Bezier method
        """
        print('now in get_bezier_xt')
    def load_PCA_xt(self):
        """
        go to the corresponding data paths to load the saved track length 
        distributions calculated using the PCA method
        """
        print('now in get_PCA_xt')


