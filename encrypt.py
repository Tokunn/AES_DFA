#!/usr/bin/env python3

import numpy as np
from tables import sbox_table, gal, Rcon, inv_sbox_table, inv_gal
gal = np.asarray(gal)
inv_gal = np.asarray(inv_gal)
np.set_printoptions(formatter={'int':'{:02x}'.format})

DEBUG_PRINT_ROUNDKEY_GEN = False
DEBUG_PRINT_ENCRIPT = True
DEBUG_PRINT_DECRYPT = False

def print_rkey(*arg):
    if DEBUG_PRINT_ROUNDKEY_GEN:
        print("[DEBUG] rkey ", *arg)

def print_dec(*arg):
    if DEBUG_PRINT_DECRYPT:
        print("[DEBUG] dec ", arg[0]+'\n', *arg[1:])

def print_enc(*arg):
    if DEBUG_PRINT_ENCRIPT:
        print("[DEBUG] enc ", arg[0]+'\n', *arg[1:])

def RotWord(word):
    # [1byte x 4]の1次元配列を受け取って左へ一つづつシフト
    return np.asarray([word[i-3] for i in range(4)])

def SubWord(word):
    return np.asarray([sbox_table[i] for i in word])

def make_extkey(key):
    # 11ラウンド分の鍵を生成
    w = []
    # 1ラウンド目の鍵はそのまま(縦に取り出す)
    for i in range(4):
        w.append(key[:,i])
    
    for i in range(4, 4*11):
        print_rkey("i", i)
        temp = w[i-1]
        print_rkey("temp", temp)
        # ラウンドの最初のみRotWordとSubWordしてRconとXORをとる
        if i%4==0:
            temp = RotWord(temp)
            print_rkey("After RotWord()", temp)
            temp = SubWord(temp)
            print_rkey("After SubWord()", temp)
            t_rcon = Rcon[int(i/4)]
            print_rkey("Rcon[i/4]", hex(t_rcon))
            temp[0] = temp[0] ^ t_rcon
            print_rkey("after xor with Rcon", temp)
        # 前のラウンド(i-4)とXOR
        t_wi4 = w[i-4]
        print_rkey("w[i-4]", t_wi4)
        wi = t_wi4 ^ temp
        print_rkey("w[i]", wi)
        # 追加　
        w.append(wi)

    # ラウンドごとに縦に成形する
    ext_keys = np.asarray(w).reshape(11, 4, 4)
    ext_keys = np.asarray([i.T for i in ext_keys])
    print_rkey(ext_keys)

    return ext_keys


def AddRoundKey(state, key):
    # Round鍵とXOR
    print_enc("Round Key Value", key)
    return state ^ key

def SubBytes(state):
    # SBox Tableで変換
    state = state.reshape(-1)
    new_state = np.zeros_like(state)
    for i in range(state.shape[0]):
        new_state[i] = sbox_table[state[i]]
    r = new_state.reshape(4,4)
    print_enc("After SubBytes", r)
    return r

def invSubBytes(state):
    # invSBox Tableで変換
    state = state.reshape(-1)
    new_state = np.zeros_like(state)
    for i in range(state.shape[0]):
        new_state[i] = inv_sbox_table[state[i]]
    r = new_state.reshape(4,4)
    print_dec("After invSubBytes", r)
    return r

def ShiftRows(state):
    # 行数に応じて横方向にシフト
    new_state = np.zeros_like(state)
    new_state[0] = state[0]
    for i in range(1, state.shape[0]):
        new_state[i, -i:] = state[i, :i]
        new_state[i, :-i] = state[i, i:]
    print_enc("After ShiftRows", new_state)
    return new_state

def invShiftRows(state):
    # 行数に応じて横方向にシフト
    new_state = np.zeros_like(state)
    new_state[0] = state[0]
    for i in range(1, state.shape[0]):
        new_state[i, :i] = state[i, -i:]
        new_state[i, i:] = state[i, :-i]
    print_enc("After invShiftRows", new_state)
    return new_state

def gmult(a, b):
    # aの最下位ビットが立ってるかを確認して、立っている場合には結果にbをかける(XOR)
    # 次のビットを見るためにBは右シフト
    # xをかけた値を作るためにAは左シフト
    # x**8 + x**4 + x**3 + x + 1 を法とするGF(2)における積
    # 8bitを超えたらオーバーフロー処理
    # x**8 = -x**4 - x**3 - x -1 を代入して計算

    c = 0
    for i in range(8):
        if (b & 0b1): # 最下位ビットが立っているなら
            c ^= a 
        b >>= 1
        a <<= 1
        if (a & 0x100): # オーバーフローしているなら
            a ^= 0b100011011 # x**8 + x**4 + x**3 + x + 1  

    return c

def xor_dot(A, B):
    # A = [[0, 0, 0, 0],
    #      [1, 1, 1, 1],
    #      [2, 2, 2, 2],
    #      [3, 3, 3, 3]]
    # B = [[0],
    #      [1],
    #      [2],
    #      [3]]
    # np.dot(A, B)の和じゃなくてXOR ver

    # Bが１行の場合には縦に直す
    if len(B.shape)==1:
        B = B.reshape(B.shape[0], 1)

    assert A.shape[1] == B.shape[0]

    r = np.zeros((A.shape[0], B.shape[1]), dtype=np.int)
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            for k in range(A.shape[1]):
                r[i,j] ^= gmult(A[i][k], B[:,j][k])
    
    # rが１列の場合には横に直す
    if r.shape[1]==1:
        r = r.reshape(-1)
    return r

def MixColumns(state):
    # 縦方向の列ごとに行列演算(足し算ではなくXOR)
    new_state = np.zeros_like(state)
    for i in range(state.shape[1]):
        new_state[:,i] = xor_dot(gal, state[:,i])
    print_enc("After MixColumns", new_state)
    return new_state
    
def invMixColumns(state):
    # 縦方向の列ごとに行列演算(足し算ではなくXOR)
    new_state = np.zeros_like(state)
    for i in range(state.shape[1]):
        new_state[:,i] = xor_dot(inv_gal, state[:,i])
    print_enc("After invMixColumns", new_state)
    return new_state

def encrypt(p_text, o_key):

    # Cast text and key
    # 縦に並べる
    state = np.asarray(p_text).reshape(4, 4).T
    # 縦に並べる
    key = np.asarray(o_key).reshape(4, 4).T
    
    # Extend Key
    ext_key = make_extkey(key)

    print_enc('-'*20+"\nStart of Round", state)
    state = AddRoundKey(state, ext_key[0])

    # 1 to 9 round
    for i in range(1, 10): # 9times
        print_enc('-'*20+"\nStart of {} Round".format(i), state)
        state = SubBytes(state)
        state = ShiftRows(state)
        state = MixColumns(state)
        state = AddRoundKey(state, ext_key[i])

    # 10 round (No MixColumns)
    print_enc('-'*20+"\nStart of {} Round".format(10), state)
    state = SubBytes(state)
    state = ShiftRows(state)
    state = AddRoundKey(state, ext_key[10])

    print('='*20+"\nENC OUTPUT\n", state)
    return state.T.reshape(-1)

def decrypt(c_text, o_key):

    # Cast text and key
    # 縦に並べる
    state = np.asarray(c_text).reshape(4, 4).T
    # 縦に並べる
    key = np.asarray(o_key).reshape(4, 4).T
    
    # Extend Key
    ext_key = make_extkey(key)

    print_dec('-'*20+"\nStart of Round", state)
    state = AddRoundKey(state, ext_key[10])
    state = invShiftRows(state)
    state = invSubBytes(state)

    for i in range(9, 0, -1): # 9times [9,8,...2,1]
        print_dec('-'*20+"\nStart of {} Round".format(i), state)
        state = AddRoundKey(state, ext_key[i])
        state = invMixColumns(state)
        state = invShiftRows(state)
        state = invSubBytes(state)

    state = AddRoundKey(state, ext_key[0])

    print('='*20+"\nDEC OUTPUT\n", state)
    return state.T.reshape(-1)

if '__main__' == __name__:
    # plane_text = 0x3243f6a8885a308d313198a2e0370734
    # original_key = 0x2b7e151628aed2a6abf7158809cf4f3c
    original_key = [0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c]

    plane_text = [0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34]
    cipher_text = encrypt(plane_text, original_key)

    # cipher_text = [0x39, 0x25, 0x84, 0x1d, 0x02, 0xdc, 0x09, 0xfb, 0xdc, 0x11, 0x85, 0x97, 0x19, 0x6a, 0x0b, 0x32]
    print(decrypt(cipher_text, original_key))

    # test = np.asarray([[0, 1, 2, 3],
    #                    [4, 5, 6, 7],
    #                    [8, 9, 0, 1],
    #                    [2, 3, 4, 5]])
    # test2 = np.asarray([1, 2, 3, 4]).reshape(4, 1)
    # print(test2)
    # print(np.dot(test, test2))
    # print(test @ test2)
    # print(xor_dot(test, test2))

    # print("ShiftRows")
    # print(ShiftRows(test))
    # print("MixColumns")
    # print(MixColumns(test))
    # print("SubBytes")
    # print(SubBytes(test))

