# Script to find regions of transposable element (TE) 'hotspots' and 'coldspots' from a BED file containing TE annotations.

# The input files are sorted BED files and the program can be called up using the following syntax:
# python3 hotspots.py INPUT_FILE GENUS SPECIES GENOME_SIZE INTERVAL

## Steps

# Calculate the mean TE count based on the manually inputted genome size.
# Iterate through each scaffold and count the number of TEs in every X base stretch (X = interval size).
# During the iteration calculate the fold change for this interval by comparing it to the baseline value (i.e. how many more TEs are in this interval relative to the baseline value?). Store the fold change, scaffold name, start and end co-ordinates as a list.
# Zip the lists together to create a nested list.
# Print the nested list into two tab-delimited files; each individual list in the nested list is printed as a separate column.

# Output should look like this:
# Scaffold Name | Start | End | Fold Change | Species

# Import the necessary modules and open the file.
import sys
import csv

# Store the command line arguments as variables.
bedFile = csv.reader(open(sys.argv[1]), delimiter='\t')
speciesName = sys.argv[2] + sys.argv[3]
genomeSize = int(sys.argv[4])

# If no interval is given in the command line, set it at the default value of 10,000bp.
if len(sys.argv) <= 5:
    interval = 10000
else:
    interval = int(sys.argv[5])

# Store the scaffold of the first line of the file as the 'current scaffold'.
currentScaffold = next(bedFile)[0]

# Count the number of TEs in the bed file. This should be equivalent to the number of lines/rows as each entry represents a separate transposon.
# TEcount needs to start at 1 as the currentScaffold = next(bedFile)[0] line means that the following for loop begins at line 2 and not line 1.
TEcount = 1
for row in bedFile:
    TEcount += 1

# Open the file again. This resets the working line to line 1, allowing the whole file to be iterated again from the start.
bedFile = csv.reader(open(sys.argv[1]), delimiter='\t')

# Set the necessary variables.
# iterCount starts at 1 as this determines what interval stretch the program is on. E.g. the first interval the user wants to look at is 0 - 10,000bp (1*iterCount), not 0 - 0 (0*iterCount).
iterCount = 1
TEcountHolder = 0
lineCount = 0
foldChangeList = []
scaffoldList = []
startList = []
endList = []
    
# Calculate the baseline by dividing TEcount by the genome size, then multiply to get the number of TEs per interval stretch.
baseline = (TEcount/genomeSize)*interval

# FOR LOOP:
# If the scaffold matches the current scaffold then count the number of TEs there are in every interval stretch. If not, then store the new scaffold as the current scaffold and begin the count again.
# New interval: reset TEcountHolder and look at the next interval (by increasing iterCount by 1).
# New scaffold: reset both TEcountHolder and iterCount to 1.
# Resetting goes to 1 and not 0 otherwise the line that the script is on when a new interval or scaffold is begun will not be counted.
for row in bedFile:
    if row[0] == currentScaffold:
        if int(row[2]) < iterCount*interval:
            TEcountHolder += 1
        # Once the interval limit has been surpassed calculate the fold change and append this to the list, along with scaffold, start and end co-ordinates.    
        else:
            # If TEcountHolder is greater than the baseline then calculate a positive fold change.
            if TEcountHolder > baseline:
                foldChange = TEcountHolder/baseline
                foldChangeList.append(foldChange)
                scaffoldList.append(currentScaffold)
                startList.append((iterCount-1)*interval)
                endList.append(iterCount*interval)
                iterCount += 1
                TEcountHolder = 1
            # If TEcountHolder is smaller than the baseline then calculate a negative fold change.
            else:
                foldChange = -(baseline/TEcountHolder)
                foldChangeList.append(foldChange)
                scaffoldList.append(currentScaffold)
                startList.append((iterCount-1)*interval)
                endList.append(iterCount*interval)
                iterCount += 1
                TEcountHolder = 1
    # When a new scaffold begins calculate and append all the necessary variables and store the new scaffold as the current scaffold and begin the count again.
    else:
        if TEcountHolder > baseline:
            foldChange = TEcountHolder/baseline
            foldChangeList.append(foldChange)
            scaffoldList.append(currentScaffold)
            startList.append((iterCount-1)*interval)
            endList.append(iterCount*interval)
            currentScaffold = row[0]
            TEcountHolder = 1
            iterCount = 1
        else:
            foldChange = -(baseline/TEcountHolder)
            foldChangeList.append(foldChange)
            scaffoldList.append(currentScaffold)
            startList.append((iterCount-1)*interval)
            endList.append(iterCount*interval)
            currentScaffold = row[0]
            TEcountHolder = 1
            iterCount = 1


# Another if block at the end to ensure the last line of the file is counted. The previous for loop does not work for the last line of the file as it is dependent on the next line to be different to the previous one to trigger the storing of the values into the list.
if TEcountHolder > baseline:
        foldChange = TEcountHolder/baseline
        foldChangeList.append(foldChange)
        scaffoldList.append(currentScaffold)
        startList.append((iterCount-1)*interval)
        endList.append(iterCount*interval)
else:
        foldChange = -(baseline/TEcountHolder)
        foldChangeList.append(foldChange)
        scaffoldList.append(currentScaffold)
        startList.append((iterCount-1)*interval)
        endList.append(iterCount*interval)

# Create the unfiltered tab-delimited output.
myFile = open(f"{speciesName}Hotspots", "w")
writeFile = csv.writer(myFile, delimiter="\t")

# Create a list that just contains the species name and is as long as the fold change list.
speciesList = [speciesName]*len(foldChangeList)

# Zip the lists together to effectively create a nested list.
finalList = zip(scaffoldList, startList, endList, foldChangeList, speciesList)

# Write the nested list to the file where each list in the nested list is printed as a separate line/entry/row.
writeFile.writerow(["Scaffold", "Start", "End", "Fold Change", "Species Name"])
for row in finalList:
    writeFile.writerow(row)

# Create a second file with only the larger fold change values.
# Need to re-create the finalList object as it is disappears after the above code.
myFile = open(f"{speciesName}Hotspots2", "w")
writeFile = csv.writer(myFile, delimiter="\t")

finalList = zip(scaffoldList, startList, endList, foldChangeList, speciesList)
writeFile.writerow(["Scaffold", "Start", "End", "Fold Change", "Species Name"])

# Only print out fold changes which are smaller than -3x or greater than 3x.
for row in finalList:
    if row[3] <= -3 or row[3] >= 3:
        writeFile.writerow(row)

print(f"Total number of TEs counted: {TEcount}")
print(f"The mean number of TEs per {interval} bases: {baseline}")
