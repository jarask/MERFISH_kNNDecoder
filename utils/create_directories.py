import os
import glob


# Set the path of which the directories should be created
path = 'C:/Users/jakob/Documents/Repositories/starfish_simulation_analysis/data/p50_sigma2'
dir_list = glob.glob('%s/*' % path + os.path.sep)
dir_list2 = []
for directory in dir_list:
    dir_list2.append(directory.split('\\')[1])

# Create directories
for directory in dir_list2:
    os.makedirs('data/p50_sigma2/%s' % directory, exist_ok=True)
