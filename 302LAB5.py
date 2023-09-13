import random
def gibbs_sampler(protein, k, t, N):
    motifs = [protein[i][random.randint(0, len(protein[i])-k):][:k] for i in range(t)]
    best_motifs = motifs
    for j in range(N):
        i = random.randint(0, t-1)
        motifs.pop(i)
        profile = create_profile(motifs)
        motifs.insert(i, profile_most_probable(protein[i], k, profile))
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = motifs
    return best_motifs

def create_profile(motifs):
    profile = []
    for i in range(len(motifs[0])):
        column = [motifs[j][i] for j in range(len(motifs))]
        profile.append([column.count('A')/len(column), column.count('C')/len(column), column.count('D')/len(column), column.count('E')/len(column), column.count('F')/len(column), column.count('G')/len(column), column.count('H')/len(column), column.count('I')/len(column), column.count('K')/len(column), column.count('L')/len(column), column.count('M')/len(column), column.count('N')/len(column), column.count('P')/len(column), column.count('Q')/len(column), column.count('R')/len(column), column.count('S')/len(column), column.count('T')/len(column), column.count('V')/len(column), column.count('W')/len(column), column.count('Y')/len(column)])
    return profile

def profile_most_probable(protein, k, profile):
    max_prob = -1
    most_probable = protein[:k]
    for i in range(len(protein)-k+1):
        pattern = protein[i:i+k]
        prob = 1
        for j in range(k):
            if pattern[j] == 'A':
                prob *= profile[j][0]
            elif pattern[j] == 'C':
                prob *= profile[j][1]
            elif pattern[j] == 'D':
                prob *= profile[j][2]
            elif pattern[j] == 'E':
                prob *= profile[j][3]
            elif pattern[j] == 'F':
                prob *= profile[j][4]
            elif pattern[j] == 'G':
                prob *= profile[j][5]
            elif pattern[j] == 'H':
                prob *= profile[j][6]
            elif pattern[j] == 'I':
                prob *= profile[j][7]
            elif pattern[j] == 'K':
                prob *= profile[j][8]
            elif pattern[j] == 'L':
                prob *= profile[j][9]
            elif pattern[j] == 'M':
                prob *= profile[j][10]
            elif pattern[j] == 'N':
                prob *= profile[j][11]
            elif pattern[j] == 'P':
                prob *= profile[j][12]
            elif pattern[j] == 'Q':
                prob *= profile[j][13]
            elif pattern[j] == 'R':
                prob *= profile[j][14]
            elif pattern[j] == 'S':
                prob *= profile[j][15]
            elif pattern[j] == 'T':
                prob *= profile[j][16]
            elif pattern[j] == 'V':
                prob *= profile[j][17]
            elif pattern[j] == 'W':
                prob *= profile[j][18]
            elif pattern[j] == 'Y':
                prob *= profile[j][19]
        if prob > max_prob:
            max_prob = prob
            most_probable = pattern
    return most_probable

def score_motifs(motifs):
    score = 0
    for i in range(len(motifs[0])):
        column = [motifs[j][i] for j in range(len(motifs))]
        score += len(column) - max(column.count('A'), column.count('C'), column.count('D'), column.count('E'), column.count('F'), column.count('G'), column.count('H'), column.count('I'), column.count('K'), column.count('L'), column.count('M'), column.count('N'), column.count('P'), column.count('Q'), column.count('R'), column.count('S'), column.count('T'), column.count('V'), column.count('W'), column.count('Y'))
    return score

# Example usage
protein = ['ACGAGATA',
           'GAGAACTA',
           'TAGAGACC',
	'CACTGAGA']
k = 4
t = 4
N = 1000
best_motifs = gibbs_sampler(protein, k, t, N)
print(best_motifs)