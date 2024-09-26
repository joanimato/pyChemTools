import discord
from discord.ext import commands
from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy
from io import BytesIO

BOT_TOKEN = "MTI4ODkwODc0OTAwNzg4NDI4OA.GmTmL8.Mj7-JS6lDO_tex_7PnAhB1k0R6JkdjcXtYyAPY"
CHANNEL_ID = 1288909652339331112


# Set up the bot with a command prefix
bot = commands.Bot(command_prefix='!', intents=discord.Intents.all())


@bot.event
async def on_ready():
    channel = bot.get_channel(CHANNEL_ID)
    await channel.send("Hello from SMILES bot!")

# Command to fetch molecule image from SMILES string
@bot.command(name='molimage')
async def molimage(ctx, *, smiles: str):
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


# Run the bot with your token
bot.run(BOT_TOKEN)
