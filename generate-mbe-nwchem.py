import sys
import os
from itertools import combinations
from math import factorial

# Read command-line arguments
args = sys.argv
file = args[1]
nbodies = int(args[2])
n1bodies = int(args[3])

# Class definitions for data structures


class TotalInput:
    def __init__(self, preamble, coords, postamble):
        self.pre = preamble
        self.coords = coords
        self.post = postamble

    @staticmethod
    def from_coords(inputlines : str):
        """Separates coordinates from the rest of the file 
        and returns a total input class member."""
        dataline = next(i for i, line in enumerate(inputlines) if find_text("geometry", line))
        endline = next(i for i, line in enumerate(inputlines[dataline:], start=dataline) if find_text("end", line))
        
        # print("inputlines pre: ", inputlines[:dataline + 1])
        return TotalInput(
            [x.strip() for x in inputlines[:dataline + 1]],
            [x.strip() for x in inputlines[dataline + 1: endline]],
            [x.strip() for x in inputlines[endline:]]
        )
    def add_to_end(self):
        """Adds an index to the end of each coordinate line
            of a TotalInput object."""
        newCoords = [f"{line.strip()}  #  {i+1}" for i, line in enumerate(self.coords)]
        # print("add to end newcoords", newCoords)
        return TotalInput(self.pre, newCoords, self.post)
    

    def change_start(self, title):
        """Changes the 'start' keyword for nwchem inputs."""
        preamb = self.pre
        startline = next(i for i, line in enumerate(preamb) if find_text("start", line))
        new_start = ["start " + title] + preamb[startline + 1:]
        return TotalInput(new_start, self.coords, self.post)
    
    def combine_input_file(self):
        """Combines fields of TotalInput into a single list of strings."""
        return self.pre + self.coords + self.post
    

    def print_list(self, filename : str):
        """Writes every line of the input data (TotalInput) into a text file."""
        with open(filename, 'w') as file:
            for line in self.combine_input_file():
                file.write(f"{line}\n")



# Functions
def find_text(text, substring):
    """Returns True if the substring is found in the text."""
    return text in substring


def combinations(n, r):
    """Calculates combinations of n objects taken r at a time."""
    return factorial(n) // (factorial(r) * factorial(n - r))

def n_list_select(n, n1, inputs):
    """Select up to n-body files."""
    if n > n1:
        raise ValueError("Impossible n-body selection")
    
    def aux(n, n1):
        return sum(combinations(n1, x) for x in range(1, n + 1))
    
    return inputs[:int(aux(n, n1))]

def sorted_powerset(s):
    """Creates powerset of any list and sorts it."""
    powerset = [[]]
    for elem in s:
        powerset += [subset + [elem] for subset in powerset]
    return sorted(powerset, key=len)[1:]

def limited_powerset(n, n1, s):
    """Returns a limited powerset."""
    return n_list_select(n, n1, sorted_powerset(s))

def reverse_powerset(s):
    """Reverses the powerset."""
    return sorted_powerset(s)[::-1]

def small_powerset(s):
    """Creates powerset without the largest term."""
    powerset = sorted_powerset(s)
    return sorted(powerset[1:], key=len)[1:]



def generate_mbe_coords(coords, n, n1):
    """Generates MBE coordinates, represented as a list of Coord objecst."""
    return [x for x in limited_powerset(n, n1, coords)]

def create_file_list(mbe_coords, filename):
    """Creates a list of filenames based on MBE coordinates."""
    def atom_number(line):
        #print("atom number", line)
        return int(line.split("#")[1].strip())

    def coord_num_ext(coords, acc=""):
        if not coords:
            return acc
        #print("coord_num_ext", coords)
        return coord_num_ext(coords[1:], acc + "_" + str(atom_number(coords[0])) + "_")

    def create_extension(coords_list):
        if not coords_list:
            return []
        return [coord_num_ext(c) for c in coords_list]

    base_name = filename[:filename.index(".nw")]
    extensions = create_extension(mbe_coords)
    
    return [f"{base_name}{ext}mbe.nw" for ext in extensions]


def add_desc_to_coord(pre, post, coords):
    """Adds file description to each MBE coordinate set."""
    return [TotalInput(pre, c, post) for c in coords]


def write_to_files(filelist, inputs):
    """Writes the MBE coordinates to the corresponding files."""
    # print("write_to_files", inputs)
    # print("write_to_files", filelist)
    for filename, coords in zip(filelist, inputs):
        coords.print_list(filename)


# Main execution
if __name__ == "__main__":
    lines = open(file).readlines()
    total_input = TotalInput.from_coords(lines)
    #print("total input pre: ", total_input.pre)
    # print("total input coord: ", total_input.coords)
    #print("total input : ", total_input.post)
    newTotal_input = total_input.add_to_end()

    mbe_coords = generate_mbe_coords(newTotal_input.coords, nbodies, n1bodies)
    #print("mbe_coords: ", mbe_coords)
    mbe_files = add_desc_to_coord(total_input.pre, total_input.post, mbe_coords)
    #for f in mbe_files:
    #    print("mbe_files: ", f.coords)

    filelist = create_file_list(mbe_coords, file)

    mbe_files_corrected = [mbe_file.change_start(f) for mbe_file, f in zip(mbe_files, filelist)]

    write_to_files(filelist, mbe_files_corrected)
