# list or csv file
cwd = os.getcwd()
doc_path=os.path.abspath('C:\\Users\\chalabaj\\Documents')
doc_listdir= os.listdir('C:\\Users\\chalabaj\\Documents')

csvfile=([xx for xx in doc_listdir if xx.endswith('.csv')])
print(doc_path, csvfile)
csvfile_path= os.path.join(doc_path,csvfile[0])
print(csvfile_path)
type(csvfile_path)
