import sys
import os

# Read command-line arguments
args = sys.argv
file = args[1]

# Class definitions for data structures
class Coord:
    def __init__(self, coords):
        self.coords = coords

class Coordinates:
    def __init__(self, coordinates):
        self.coordinates = coordinates

class MBECoord:
    def __init__(self, coords):
        self.coords = coords

class TxtInput:
    def __init__(self, text):
        self.text = text

class TotalInput:
    def __init__(self, preamble, coords, postamble):
        self.preamble = preamble
        self.coords = coords
        self.postamble = postamble

# Functions
def find_text(text, substring):
    """Returns True if the substring is found in text."""
    return substring.find(text) != -1

def parse_coords(inputlines):
    """Separates coordinates from the rest of the file."""
    dataline = next(i for i, line in enumerate(inputlines) if find_text(" $DATA", line))
    endline = next(i for i, line in enumerate(inputlines[dataline:], start=dataline) if find_text(" $END", line))
    
    return TotalInput(
        TxtInput(inputlines[:dataline + 2]),
        Coord(inputlines[dataline + 3: dataline + endline - 1]),
        TxtInput(inputlines[dataline + endline:])
    )

def convert_to_coord(textcords):
    """Converts text coordinates to a float 2D array."""
    arrcoord = textcords.coords
    splitarr = [list(filter(lambda x: x != "", line.split())) for line in arrcoord]
    c = [[float(x) for x in coord[2:5]] for coord in splitarr]
    return Coordinates(c)

def parse_desc(inp_file):
    """Extracts the description part from the input file."""
    dataline = next(i for i, line in enumerate(inp_file) if find_text(" $DATA", line))
    return inp_file[:dataline]

def sorted_powerset(s):
    """Generates a sorted powerset of the input list."""
    powerset = [[]]
    for elem in s:
        powerset += [subset + [elem] for subset in powerset]
    return sorted(powerset, key=len)[1:]

def reverse_powerset(s):
    """Returns the reversed powerset."""
    return sorted_powerset(s)[::-1]

def small_powerset(s):
    """Generates a powerset without the largest term."""
    powerset = sorted_powerset(s)
    return sorted(powerset[1:], key=len)[1:]

def add_to_end(coords):
    """Adds an index to the end of each coordinate line."""
    return Coord([f"{line}  !  {i+1}" for i, line in enumerate(coords.coords)])

def generate_mbe_coords(coords):
    """Generates the MBE coordinates."""
    return MBECoord([Coord(x) for x in sorted_powerset(coords.coords)])

def create_file_list(mbe_coords, filename):
    """Creates a list of filenames based on MBE coordinates."""
    def atom_number(line):
        return int(line.split("!")[1].strip())

    def coord_num_ext(coords, acc=""):
        if not coords:
            return acc
        return coord_num_ext(coords[1:], acc + str(atom_number(coords[0])))

    def create_extension(coords_list):
        if not coords_list:
            return []
        return [coord_num_ext(c.coords) for c in coords_list]

    base_name = filename[:filename.index(".inp")]
    extensions = create_extension(mbe_coords.coords)
    
    return [f"{base_name}{ext}mbe.inp" for ext in extensions]

def create_atom_array(mbe_coords):
    """Creates an array of atom indices from MBE coordinates."""
    def atom_number(line):
        return int(line.split("!")[1].strip())

    def construct_one_array(coords, acc=[]):
        if not coords.coords:
            return acc
        return construct_one_array(Coord(coords.coords[1:]), acc + [atom_number(coords.coords[0])])

    return [construct_one_array(c) for c in mbe_coords.coords]

def add_desc_to_coord(pre, post, coords):
    """Adds file description (preamble and postamble) to each MBE coordinate set."""
    return [TotalInput(pre, c, post) for c in coords.coords]

def combine_input_file(input_data):
    """Combines the fields of TotalInput into a single list of strings."""
    return input_data.preamble.text + input_data.coords.coords + input_data.postamble.text

def print_list(filename, input_data):
    """Writes a TotalInput object to a file."""
    with open(filename, 'w') as file:
        for line in combine_input_file(input_data):
            file.write(f"{line}\n")

def write_to_files(filelist, mbe_coords):
    """Writes the MBE coordinates to the corresponding files."""
    for filename, coords in zip(filelist, mbe_coords):
        print_list(filename, coords)

def create_mbe(inpfile):
    """Main function to generate all possible many-body input files."""
    file = inpfile
    lines = open(file).readlines()

    pre, coor, post = parse_coords(lines)
    coor = add_to_end(coor)

    mbe_coords = generate_mbe_coords(coor)
    mbe_files = add_desc_to_coord(pre, post, mbe_coords)
    filelist = create_file_list(mbe_coords, file)

    write_to_files(filelist, mbe_files)

# Execute the MBE creation process
if __name__ == "__main__":
    lines = open(file).readlines()
    pre, coor, post = parse_coords(lines)
    coor = add_to_end(coor)

    mbe_coords = generate_mbe_coords(coor)
    mbe_files = add_desc_to_coord(pre, post, mbe_coords)
    filelist = create_file_list(mbe_coords, file)

    write_to_files(filelist, mbe_files)
