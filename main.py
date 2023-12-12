import galois  # GF operations
from Crypto.Random import get_random_bytes  # key and npub generation
from Crypto.Cipher import AES
from Crypto.Util.number import bytes_to_long, long_to_bytes  # for operations on GF
from Crypto.Util.strxor import strxor  # xor two byte strings
from functools import reduce


def rho(x, st):
    gfmult1 = long_to_bytes(
        (GF(3) * GF(bytes_to_long(st))).item()
    )
    gfmult2 = long_to_bytes(
        (GF(2) * GF(bytes_to_long(st))).item()
    )
    y = strxor(x, gfmult1)
    st_prim = strxor(x, gfmult2)
    return y, st_prim


def rho_inverse(y, st):
    gfmult = long_to_bytes(
        (GF(3) * GF(bytes_to_long(st))).item()
    )
    x = strxor(y, gfmult)
    st_prim = strxor(y, st)
    return x, st_prim


def delta_A(i, L, A_star_len, a):
    if i+1 == a and A_star_len < 16:
        return long_to_bytes(
            (GF(3) * GF(7) * GF(2)**(i-1) * GF(bytes_to_long(L))).item()
        )
    else:
        return long_to_bytes(
            (GF(3) * GF(2)**i * GF(bytes_to_long(L))).item()
        )


def delta_M(i, L, M_star_len, l):
    if i in (l, l+1) and M_star_len == 16:
        return long_to_bytes(
            (GF(7) * GF(2)**(i-1) * GF(bytes_to_long(L))).item()
        )
    elif i in (l, l+1) and M_star_len < 16:
        return long_to_bytes(
            (GF(7)**2 * GF(2) ** (i - 1) * GF(bytes_to_long(L))).item()
        )
    else:
        return long_to_bytes(
            (GF(2)**i * GF(bytes_to_long(L))).item()
        )


def delta_C(i, L, M_star_len, l):
    # Using the fact that delta_C = GF(3)**2 delta_M
    deltaM = delta_M(i, L, M_star_len, l)
    return long_to_bytes(
        (GF(3)**2 * GF(bytes_to_long(deltaM))).item()
    )


def gen_subkey(k):
    cipher = AES.new(k, AES.MODE_ECB)
    L = cipher.encrypt(b'\x00' * 16)  # Did author mean 128 bits of ascii zeroes of 0 bytes?
    return L


def onezero_pad(b):
    # Byte object is immutable, so we need to create a new one. We can't append to existing
    pad_len = 16 - len(b)  # Amount of bytes to pad
    if pad_len != 0:
        b += b'\x80'  # '1'
        for _ in range(pad_len - 1):  # Fill with zeroes padding to 128 bits
            b += b'\x00'
    return b


def generate_iv(k, A, param, npub, L):
    # subkey == L

    # Separate associated data into blocks

    A = [A[i:i+16] for i in range(0, len(A), 16)]
    a = len(A)  # required for masking functions. Amount of blocks
    a_star_len = len(A[-1])  # required for masking functions

    # apply '10' padding for the last block if needed
    A[-1] = onezero_pad(A[-1])

    # Step 1 - W'[0]
    cipher = AES.new(k, AES.MODE_ECB)
    W_0 = cipher.encrypt(
        strxor((npub + param), delta_A(0, L, a_star_len, a))
    )

    # Step 2 - AA[i]
    AA = []
    for i in range(a):
        AA.append(
            strxor(A[i], delta_A(i, L, a_star_len, a))
        )
    # Step 3 Z[i]
    Z = []
    for i in range(a):
        Z.append(
            cipher.encrypt(AA[i])
        )
    # Step 4 - W' for i != 0
    W_prim = []
    W_prim.append(  # W[0] in paper needs special var, since from paper W[1] is W[0] in Python
        strxor(Z[0], W_0)
    )
    for i in range(1, a):
        W_prim.append(
            strxor(Z[i], W_prim[i-1])
        )

    return W_prim[-1]


def generate_tagged_ct(k, M, IV):
    # Separate plaintext message into blocks
    M = [M[i:i + 16] for i in range(0, len(M), 16)]
    l = len(M) + 1  # required for masking functions. Amount of blocks
    M_star_len = len(M[-1])  # required for masking functions
    # Apply padding to last block and generate it based on the provided method from the paper
    M[-1] = onezero_pad(M[-1])
    M[-1] = reduce(lambda x, y: strxor(x, y), M)
    M.append(M[-1])

    # Step 1 - set W_0
    W_0 = IV

    # Step 2 - MM
    MM = []
    for i in range(l):
        MM.append(
            strxor(M[i], delta_M(i, L, M_star_len, l))
        )


if __name__ == "__main__":
    # AES block_size is usually 16?
    # For COLM0 tau = 0, l_r = 128
    # galois.GF creates FieldArray subclass which is a subclass of ndarray from numpy. Debugger can't comprehend that...
    GF = galois.GF(2 ** 128, irreducible_poly="x^128 + x^7 + x^2 + x + 1")

    # Generate params
    k = get_random_bytes(16)
    npub = get_random_bytes(8)
    #  tau = 0, next 8 bits represent l_r = 128, then 40 bits of zeroes
    param = b'\x00' * 2 + b'\x80' + b'\x00' * 5
    associated_data = b"Associated data is visible"
    plaintext = b"This should be kept secret"

    L = gen_subkey(k)
    IV = generate_iv(k, associated_data, param, npub, L)
    C = generate_tagged_ct(k, plaintext, IV)
    C = b''.join(C)
    pass

