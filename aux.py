nextdir = '/Users/jjgomezcadenas/NEXT'
icdir = os.path.join(nextdir,"IC") 
icbindir = os.path.join(icdir,"bin") 

os.environ['ICTDIR'] = icdir 
os.environ["PATH"] = os.environ["PATH"] + ':' + icbindir

print(os.environ["PATH"])
