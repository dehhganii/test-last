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
    return text

"""
text = gen_solid_string(150, 4)
text = gen_indet_string(text, 4, 20)
#print(text)
pattern = gen_solid_string(4, 4)
pattern = gen_indet_string(pattern, 4, 2)
#print(pattern)
"""
def generate_text_pattern(text_length, text_indet_letters, pattern_length, pattern_indet_letters, alphabet):
    text = gen_solid_string(text_length, alphabet)
    text = gen_indet_string(text, alphabet, text_indet_letters)
    pattern = gen_solid_string(pattern_length, alphabet)
    pattern = gen_indet_string(pattern, alphabet, pattern_indet_letters)
    f = open("/Users/hossein/Documents/Final_result/demofile_pattern", "w")
    for i in pattern:
        f.write(str(i) + ",")
    f.close()

    f2 = open("/Users/hossein/Documents/Final_result/demofile_text", "w")
    for i in text:
        f2.write(str(i) + ",")
    f2.close()

"""
for i in range(1, 300):
    generate_text_pattern(500 * i, 50, 20, 5, 4)
    subprocess.call("./a.out")
"""


print(gen_power_set(gen_prime(4)))