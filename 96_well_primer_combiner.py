#! /usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import math
import string


def run(args):
    filename = args.input # these match the "dest": dest="input"
    output_filename = args.output # from dest="output"


    # Do stuff
    # Read in data
    df = pd.read_csv(filename, sep=',',names=('plateN', 'V2', 'V3', 'V4', 'V5', 'wellPos', 'seqName', 'V8',
                       		 'V9', 'V10', 'V11',
                       		 'V12', 'V13', 'V14',
           			 'V15',
                       		 'V16', 'V17', 'V18',
                       		 'V19'), quotechar='"')


    # Sort the data follow the plate order
    df['order'] = [0]*len(df)
    order = []
    totalRow = df.shape[0]
    for i in range(0, totalRow):
        order.append(i)
    df['order'] = order

    sort_by_plateN_order = df.sort_values(['plateN', 'order'])

        	        #print(sort_by_plateN_order)

    # Seperate the whole dataframe and build the dataset for each plate as name[0:8]
    numOfPlate = math.floor(totalRow / 96)  # Default will be an even number of plates
        	        #print(numOfPlate)
    name = [''] * numOfPlate
    for i in range(0, numOfPlate):
        name[i] = sort_by_plateN_order.iloc[i*96:(i*96+96), :]
                #   print(name[i][['plateN', 'wellPos', 'seqName']])



    # Build a (numOfPlate//2)*97 2D array new table for (numOfPlate//2) new plate to save their data.
    # Each plate has 97 cells for a new plate ID and data
    newPlate = [[] * 97 for i in range(numOfPlate//2)]
    for l in range(0, numOfPlate // 2):
        newPlate[l].append('plate' + str(2 * l + 1) + '_' + str(2 * l + 2))

    # Seperate all the number of plates into pairs, and combine every pair of plates into a new plate
    for i in range(0, numOfPlate, 2):
        i = i
        j = i+1
        indexNewPlate = i//2
        vector1 = []
        vector2 = []
        # Get the first plate of the pair into a 8*6 matrix after mix the forward and reverse primers together
        for a in range(0, 96, 2):
            if name[i]['seqName'].values[a] == 'Null':
                vector1.append('Null')
            else:
                vector1.append(name[i]['seqName'].values[a][:-2])
        # Get the second plate of the pair into a 8*6 matrix after mix the forward and reverse primers together
        for b in range(0, 96, 2):
            if name[j]['seqName'].values[b] == 'Null':
                vector2.append('Null')
            else:
                vector2.append(name[j]['seqName'].values[b][:-2])

        matrix1 = np.array(vector1).reshape(8, 6)

        matrix2 = np.array(vector2).reshape(8, 6)
   	        #matrix2 = np.array(vector2).reshape(8, 6)
    		# Horizontally combine the two matrix of plates together and form a 8*12 matrix
        newMatrix = np.hstack([matrix1, matrix2])
        # Collect the data in the matrix column by column into one array
        for m in range(0, 12):
            for n in range(0, 8):
                newPlate[indexNewPlate].append(newMatrix.item((n, m)))

	#print(newPlate[3])

    col1 = []
    col1.append('Well')
    for m in range(0, 12):
        for n in range(0, 8):
            letter = chr(n+65)
            wellNum = str(m+1)
            col1.append('1'+letter+wellNum.zfill(2))


    table = pd.DataFrame.from_dict(col1)
    for i in range(0, numOfPlate//2):
        table[str(i+1)] = newPlate[i]
    #print(table)

    table.to_csv(output_filename, header=None, index=None,
             sep=',', mode='a')




def main():
    parser=argparse.ArgumentParser(description="Combine the pairwise forward and reverse primers from two wells of 96 well plate into one, thus merge two plates into one")
    parser.add_argument("-in",help="input .csv file with columns named 'plateN', 'wellPos', 'seqName'; default number of plates is even and plates start from plate 1;default that data starts from the first line (doesn't have column name)" ,dest="input", type=str, required=True)
    parser.add_argument("-out",help="merged plates output .csv filename" ,dest="output", type=str, required=True)
    #parser.add_argument("-q",help="Quality score to fill in (since fasta doesn't have quality scores but fastq needs them. Default=I" ,dest="quality_score", type=str, default="I")
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
