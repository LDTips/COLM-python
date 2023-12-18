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

def remove_onezero_pad(b):
    pad_len = len(b)
    while pad_len > 0:
        if b.find(b'\x80'+b'\x00'* (pad_len - 1)) > 0:
            b = b[:-pad_len]
            return b
        pad_len-=1
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


def generate_tagged_ct(k, M, IV, L):
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
            strxor(M[i], delta_M(i+1, L, M_star_len, l-1))
        )

    # Step 3 - X
    X = []
    cipher = AES.new(k, AES.MODE_ECB)
    for i in range(l):
        X.append(
            cipher.encrypt(MM[i])
        )

    # Step 4 - Y and W
    Y, W = [], []
    x, st = rho(X[0], W_0)
    Y.append(x)
    W.append(st)

    for i in range(1, l):
        x, st = rho(X[i], W[i-1])
        Y.append(x)
        W.append(st)

    # Step 5 - CC
    CC = []
    for i in range(l):
        CC.append(
            cipher.encrypt(Y[i])
        )

    # Step 6 - C (ciphertext)
    C = []
    for i in range(l):
        C.append(
            strxor(CC[i], delta_C(i+1, L, M_star_len, l-1))
        )
    C[l-1] = C[l-1][:M_star_len]
    return C

def decrypt_tagged_ct(k, C, IV, L):
    # Separate ciphertext into blocks
    C = [C[i:i + 16] for i in range(0, len(C), 16)]
    l = len(C) - 1
    C_last_len = len(C[-1])

    # Step 1 - set W_0
    W_0 = IV

    # Step 2 - CC
    CC = []
    for i in range(l):
        CC.append(
            strxor(C[i], delta_C(i+1, L, C_last_len, l))
        )

    # Step 3 - Y
    Y = []
    cipher = AES.new(k, AES.MODE_ECB)
    for i in range(l):
        Y.append(
            cipher.decrypt(CC[i])
        )

    # Step 4 - X and W
    X, W = [], []
    x, st = rho_inverse(Y[0], W_0)
    X.append(x)
    W.append(st)
    for i in range(1, l):
        x, st = rho_inverse(Y[i], W[i-1])
        X.append(x)
        W.append(st)

    # Step 5 - MM
    MM = []
    for i in range(l):
        MM.append(
            cipher.decrypt(X[i])
        )

    # Step 6 - M (message)
    M = []
    for i in range(l):
        M.append(
            strxor(MM[i], delta_M(i+1, L, C_last_len, l))
        )

    # Step 7 - Add M[l+1]
    M_star = reduce(lambda x, y: strxor(x, y), M)
    M_last = M[-1]
    M[-1] = remove_onezero_pad(M_star)

    return M, M_last, W[-1]

def verify_message(k, M, M_last, W_last, C_last, L):
    M_last_len = len(M[-1])
    l = len(M) + 1
    # Step 1 - MM
    MM = strxor(M_last, delta_M(l, L, M_last_len, l))

    # Step 2 - X
    cipher = AES.new(k, AES.MODE_ECB)
    X = cipher.encrypt(MM)

    # Step 3 - Y and W
    Y, W = rho(X, W_last)

    # Step 4 - CC
    CC = cipher.encrypt(Y)

    # Step 5 - C_prim
    C_prim = strxor(CC, delta_C(l, L, M_last_len, l))
    C_prim = C_prim[:M_last_len]

    if C_prim == C_last:
        return "Verification successful"
    else:
        return "Verification failed"


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
    C = generate_tagged_ct(k, plaintext, IV, L)
    C_last = C[-1]
    ciphertext = b''.join(C)
    print(ciphertext)

    M, M_last, W_last = decrypt_tagged_ct(k, ciphertext, IV, L)
    message = b''.join(M)
    print(message)

    print(verify_message(k, M, M_last, W_last, C_last, L))

    pass

