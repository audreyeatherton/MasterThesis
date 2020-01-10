# Python code to make things

def picardAddReadGroups(runname):
  # Load libraries
  import os
  
  # Make loop
  for file in os.listdir():
    # Import data from bam files
    path = os.path.abspath(file)
    filename = file

    # Build necessary things
    lanenumber = file[file.find("L") + 1]
      # Will have 1 for yoda10 and 9 for yoda9
      # If working with double digits, change to following
    # lanenumber = file[file.find("L") + 1:file.find("L") + 2]
    runnumber = path[path.find(runname) + len(runname)]

    # Name variables for picard
    i = path
    o = "../picardSortRG/" + file
    rgid = lanenumber + runnumber
    rglb = "lib1"
    rgpl = "illumina"
    rgpu = runname
    rgsm = file[0:4]
	
    # Run picard
    pic = ("picard-tools AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (i, o, rgid, rglb, rgpl, rgpu, rgsm))
    os.system(pic)
    
picardAddReadGroups("yoda")

