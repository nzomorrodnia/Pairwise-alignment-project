def compare(input1, input2):
  '''Compares two inputs, outputs list of lists'''
  alphabet = ['A','C','G','T']
  sequences = [input1,input2]
  
  # Filtering inputs
  for seq in sequences:
    if type(seq) is not str:
       return "Please input two strings for alignment"
    seq = seq.upper()
    for n in seq:
      if n not in alphabet:
        return "Please input a valid sequence"
  
  # Calculating length and building matrix
  alignment_matrix = []
  for i in sequences[1]:
    row = []
    for j in sequences[0]:
      elem = 0
      if i==j:
        elem+=1
      row.append(elem)
    alignment_matrix.append(row)
  
  return alignment_matrix

  
for row in (compare('AGCCA','AGGACT')):
  print(row)





