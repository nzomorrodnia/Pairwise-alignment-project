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
  i = len(input1)
  j = len(input2)
  

  
print(compare('AG1CCA','AGGACT'))





