import random
from itertools import chain, combinations
import subprocess

def gen_prime(alphabet_size):
    primes = []
    count = 0
    for num in range(2, 100):
        if num > 1:
            for i in range(2, num):
                if (num % i) == 0:
                    break
            else:
                primes.append(num)
                count += 1
                if count == alphabet_size:
                    return list(primes)


def gen_power_set(primes):
    s = list(primes)
    multiply = []
    power_set = list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))
    for i in range(2 ** len(primes)):
        if len(power_set[i]) > 1:
            count = 1
            for j in range(len(power_set[i])):
                count *= power_set[i][j]
            multiply.append(count)
    return list(multiply)


def gen_solid_string(length, alphabet_size):
    text = []
    alphabet = gen_prime(alphabet_size)
    for i in range(length):
        text.append(alphabet[random.randint(0, alphabet_size - 1)])
    return text


def gen_indet_string(text, alphabet_size, num):
    alphabet = gen_prime(alphabet_size)
    indet_list = gen_power_set(alphabet)
    for i in range(0, num):
        pos = random.randint(0, len(text) - 1)
        if text[pos] in alphabet:
            text[pos] = indet_list[random.randint(0, len(indet_list) - 1)]
        else:
            pos = random.randint(0, len(text) - 1)
            text[pos] = indet_list[random.randint(0, len(indet_list) - 1)]
        #text[len(text) - num + i] = indet_list[random.randint(0, len(indet_list) - 1)]
    return text




def generate_text_pattern(text_length, text_indet_letters, pattern_length, pattern_indet_letters, alphabet):

    text = gen_solid_string(text_length, alphabet)
    text = gen_indet_string(text, alphabet, text_indet_letters)
    pattern = gen_solid_string(pattern_length, alphabet)
    pattern = gen_indet_string(pattern, alphabet, pattern_indet_letters)

    """
    text = []
    pattern = []
    for j in range(1, 500 * i):
        if j % 200 == 0:
            text.append(35)
        else:
            text.append(2)

    for j in range(20* i):
        pattern.append(2)
    
    """
    #print(text)
    #print(pattern)

    f = open("/demofile_pattern", "w")
    for i in pattern:
        f.write(str(i) + ",")
    f.close()

    f2 = open("/demofile_text", "w")
    for j in text:
        f2.write(str(j) + ",")
    f2.close()

"""
for i in range(1, 100):
    generate_text_pattern(i, 0, 20, 0, 8)
    subprocess.call("./a.out")
"""
def new_text_and_pattern(n, m):
    patt = [2] * m
    tex = patt[1:] + [15]
    text = tex * n
    #print(patt)
    #print(text)
    f = open("/demofile_pattern", "w")
    for i in patt:
        f.write(str(i) + ",")
    f.close()

    f2 = open("/demofile_text", "w")
    for j in text:
        f2.write(str(j) + ",")
    f2.close()

"""
for i in range(10, 1900, 10):
    new_text_and_pattern(i, 5 + int(i/10))
    subprocess.call(("./a.out"))
"""
#for i in range(1, 300):
 #   generate_text_pattern(500*i,int(i/50), 5+int(i/2), 0, 4)
  #  subprocess.call("./a.out")

"""
text = []
pattern = []
for i in range(1, 100000):
    if i % 200 == 0:
        text.append(35)
    else:
        text.append(2)

for i in range(200):
    pattern.append(2)

f = open("/Users/hossein/Documents/generate_test_cases/demofile_pattern", "w")
for i in pattern:
    f.write(str(i) + ",")
f.close()

f2 = open("/Users/hossein/Documents/generate_test_cases/demofile_text", "w")
for i in text:
    f2.write(str(i) + ",")
f2.close()

subprocess.call("./a.out")
"""
#print(sorted(gen_power_set(gen_prime(4))))

for i in range(1, 200):
    n = 500 * i
    m = 5 + int(0.3 * i)
    k1 = int(i / 30)
    k2 = 3
    generate_text_pattern(n, k1, m, k2, 4)
    subprocess.call("./a.out")


