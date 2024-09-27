import discord
from discord.ext import commands
from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy
from io import BytesIO
import requests 

BOT_TOKEN = "MTI4ODkwODc0OTAwNzg4NDI4OA.G1qEOW.noWhQlb7i6qsIDrX7zJ2NcA7IIoVwjqcivHlPA"
CHANNEL_ID = 1288909652339331112


# Set up the bot with a command prefix
bot = commands.Bot(command_prefix='!', intents=discord.Intents.all())


@bot.event
async def on_ready():
    channel = bot.get_channel(CHANNEL_ID)
    await channel.send("Hello from SMILES bot!")

# Command to fetch molecule image from SMILES string
@bot.command(name='molimage')
async def molimage(ctx, smiles: str):
    try:
        # Convert SMILES string to a RDKit molecule object
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            await ctx.send("Invalid SMILES string.")
            return
        
        # Generate a PNG image of the molecule
        img = Draw.MolToImage(mol)

        # Save the image to a BytesIO object
        with BytesIO() as image_binary:
            img.save(image_binary, 'PNG')
            image_binary.seek(0)
            # Send the image in the Discord channel
            await ctx.send(file=discord.File(fp=image_binary, filename='molecule.png'))

    except Exception as e:
        await ctx.send(f"An error occurred: {e}")

@bot.command(name = "molname")
async def molname(ctx, *, smiles: str): 
    try: 
        compound = pubchempy.get_compounds(smiles, namespace="smiles")
        await ctx.send(f"IUPAC name: {compound[0].iupac_name}")
    except Exception as e:
        await ctx.send(f"An error occured: {e}")


@bot.command(name = "nameToSmiles")
async def nameToSmiles(ctx, name):
    """
    Converts a common or IUPAC name to an RDKit Mol object
    and return the smiles string.
    
    Parameters:
    - name (str): The common or IUPAC name of the molecule.
    
    Returns:
    - SMILES string corresponding to the molecule, or None if not found.
    """
    try:
        # PubChem API to get compound information using the name
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for bad responses

        # Extract the SMILES from the response
        data = response.json()
        smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
        
        # Convert the SMILES to an RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            await ctx.send(f"Unable to parse SMILES for molecule: {name}")
            return None
        
        await ctx.send(smiles) 
        return mol

    except requests.exceptions.HTTPError as http_err:
         await ctx.send(f"HTTP error occurred: {http_err}")
    except Exception as err:
         await ctx.send(f"Error occurred: {err}")
    
    return None


@bot.command(name = "nameToImage")
async def nameToImage(ctx, name):
    """
    Converts a common or IUPAC name to an RDKit Mol object
    and return the image.
    
    Parameters:
    - name (str): The common or IUPAC name of the molecule.

    Returns:
    - mol object. outputs image to discord. 
    """
    try:
        # PubChem API to get compound information using the name
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for bad responses

        # Extract the SMILES from the response
        data = response.json()
        smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
        
        # Convert the SMILES to an RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            await ctx.send(f"Unable to parse SMILES for molecule: {name}")
            return None
        
        await molimage(ctx, smiles)
        return mol

    except requests.exceptions.HTTPError as http_err:
         await ctx.send(f"HTTP error occurred: {http_err}")
    except Exception as err:
         await ctx.send(f"Error occurred: {err}")
    
    return None

# Run the bot with your token
bot.run(BOT_TOKEN)
