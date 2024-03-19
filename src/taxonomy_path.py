
# Name: taxonomy_path.py
# Description: This code contains a class definition that
#              will allow us to compare taxonomy paths.

class TaxPath:
    INPUT_SOURCES = ['silva', 'ebi']

    def __init__(self, path, input_source, silva_dict=None):

        if input_source not in self.INPUT_SOURCES:
            raise ValueError(f"input_source must be one of {self.INPUT_SOURCES}")
        self.input_source = input_source
        self.full_path = path

        self.domain_ = ""
        self.phylum_ = ""
        self.class_ = ""
        self.order_ = ""
        self.family_ = ""
        self.genus_ = ""
        self.rank_level = ""

        # fill in each of the levels based on the input source
        if self.input_source == "silva":
            assert silva_dict is not None, "silva_dict must be provided for silva input_source"
            self.parse_silva_taxonomy(self.full_path, silva_dict)
        elif self.input_source == "ebi":
            self.parse_ebi_taxonomy(self.full_path)
    
    def parse_silva_taxonomy(self, full_path, silva_dict):
        """ parse the SILVA taxonomy and identiy each level """

        # split the full path into individual levels and check what rank they are
        lineage_split = full_path.split(";")
        individual_levels = [";".join(lineage_split[:i+1]) + ";" for i in range(len(lineage_split)-1)]

        for i, level in enumerate(individual_levels):
            assert level in silva_dict, f"level {level} not found in silva_dict"
            curr_rank = silva_dict[level]

            if curr_rank == "domain":
                self.domain_ = lineage_split[i]
                self.rank_level = "domain"
            elif curr_rank == "phylum":
                self.phylum_ = lineage_split[i]
                self.rank_level = "phylum"
            elif curr_rank == "class":
                self.class_ = lineage_split[i]
                self.rank_level = "class"
            elif curr_rank == "order":
                self.order_ = lineage_split[i]
                self.rank_level = "order"
            elif curr_rank == "family":
                self.family_ = lineage_split[i]
                self.rank_level = "family"
            elif curr_rank == "genus":
                self.genus_ = lineage_split[i]
                self.rank_level = "genus"

        assert self.rank_level != "", "rank_level was not set"

    def parse_ebi_taxonomy(self, full_path):
        """ 
        parse the EBI taxonomy classification and identiy each level 
        .e.g. sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria;o__Micrococcales;f__Brevibacteriaceae;g__Brevibacterium
        """
            
        # split the full path into individual levels and check what rank they are
        lineage_split = full_path.split(";")
        for level in lineage_split:
            if "sk__" in level and len(level) > 3:
                self.domain_ = level.replace("sk__", "").replace("_",  " ")
                self.rank_level = "domain"
            elif "p__" in level and len(level) > 3:
                self.phylum_ = level.replace("p__", "").replace("_",  " ")
                self.rank_level = "phylum"
            elif "c__" in level and len(level) > 3:
                self.class_ = level.replace("c__", "").replace("_",  " ")
                self.rank_level = "class"
            elif "o__" in level and len(level) > 3:
                self.order_ = level.replace("o__", "").replace("_",  " ")
                self.rank_level = "order"
            elif "f__" in level and len(level) > 3:
                self.family_ = level.replace("f__", "").replace("_",  " ")
                self.rank_level = "family"
            elif "g__" in level and len(level) > 3:
                self.genus_ = level.replace("g__", "").replace("_",  " ")
                self.rank_level = "genus"

                # handle special case for E. coli
                if self.genus_ == "Escherichia":
                    self.genus_ = "Escherichia-Shigella"

        #assert self.rank_level == "genus", "rank_level was not set"

    def get_fingerprint(self):
        """ 
        return: fingerprint of the taxonomic rank 
        
        fingerprint - concatenation of order, family, genus 
        e.g. Nitrosocaldales;Nitrosocaldaceae;uncultured;
        """
        # return f"{self.order_};{self.family_};{self.genus_};"
        return f"{self.genus_};"

    def get_fullpath(self):
        """ return the full path of the taxonomic rank """
        if self.input_source == "silva":
            return self.full_path
        elif self.input_source == "ebi":
            return f"{self.domain_};{self.phylum_};{self.class_};{self.order_};{self.family_};{self.genus_};"
        raise ValueError("input_source must be one of ['silva', 'ebi']")