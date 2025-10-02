#!/usr/bin/env python3
# -*- coding: utf-8 -*-

EPS = 1e-12
PRINT_DECIMALS = 6

# ---------------------------
# Pretty printing
# ---------------------------
def printMatrix(mx):
    n = len(mx)
    print("𐌲" + "        " * n + "ㄱ")
    for i in range(n):
        print("|", end=" ")
        for j in range(n):
            print(f"{mx[i][j]:.{PRINT_DECIMALS}f}", end=" ")
        print("|")
    print("ㄴ" + "        " * n + "┘")

# ---------------------------
# Basic helpers
# ---------------------------
def getTransposeMatrix(matrix):
    n = len(matrix)
    return [[matrix[j][i] for j in range(n)] for i in range(n)]

def getMinorMatrix(matrix, r, c):
    return [row[:c] + row[c+1:] for row in (matrix[:r] + matrix[r+1:])]

def getMatrixDeterminant(matrix):
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    if n == 2:
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]

    det = 0.0
    for c in range(n):
        sign = -1.0 if (c & 1) else 1.0
        det += sign * matrix[0][c] * getMatrixDeterminant(getMinorMatrix(matrix, 0, c))
    return det

def has_inverse(matrix):
    det = getMatrixDeterminant(matrix)
    return abs(det) >= EPS

# ---------------------------
# Inverse (Determinant/Adjugate)
# ---------------------------
def inverse_by_determinant(m):
    n = len(m)
    det = getMatrixDeterminant(m)
    if abs(det) < EPS:
        raise ValueError("행렬식이 0입니다. 역행렬이 존재하지 않습니다.")

    if n == 1:
        return [[1.0 / m[0][0]]]

    if n == 2:
        return [[ m[1][1]/det, -m[0][1]/det],
                [-m[1][0]/det,  m[0][0]/det]]

    cof = []
    for r in range(n):
        row = []
        for c in range(n):
            minor = getMinorMatrix(m, r, c)
            sign = -1.0 if ((r + c) & 1) else 1.0
            row.append(sign * getMatrixDeterminant(minor))
        cof.append(row)

    adj = getTransposeMatrix(cof)

    for i in range(n):
        for j in range(n):
            adj[i][j] /= det

    return adj

# ---------------------------
# Inverse (Gauss–Jordan)
# ---------------------------
def inverse_by_gauss_jordan(a):
    n = len(a)
    aug = [a[i][:] + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    for col in range(n):
        pivot_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot_row][col]) < EPS:
            raise ValueError("소거 과정 중 특이 행렬 발견 (피벗이 0).")

        if pivot_row != col:
            aug[col], aug[pivot_row] = aug[pivot_row], aug[col]

        pivot = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= pivot

        for r in range(n):
            if r == col:
                continue
            factor = aug[r][col]
            if abs(factor) > EPS:
                for j in range(2 * n):
                    aug[r][j] -= factor * aug[col][j]

    return [row[n:] for row in aug]

# ---------------------------
# Checks
# ---------------------------
def are_same(m1, m2, tol=1e-6):
    n = len(m1)
    for i in range(n):
        for j in range(n):
            if abs(m1[i][j] - m2[i][j]) > tol:
                return False
    return True

def multiply(A, B):
    n = len(A)
    m = len(B[0])
    p = len(B)
    assert len(A[0]) == p
    C = [[0.0]*m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            s = 0.0
            for k in range(p):
                s += A[i][k] * B[k][j]
            C[i][j] = s
    return C

def is_identity(M, tol=1e-6):
    n = len(M)
    for i in range(n):
        for j in range(n):
            expected = 1.0 if i == j else 0.0
            if abs(M[i][j] - expected) > tol:
                return False
    return True

# ---------------------------
# I/O
# ---------------------------
def read_square_matrix():
    k_raw = input("정방행렬의 차수를 입력하세요: ").strip()
    try:
        k = int(k_raw)
    except ValueError:
        raise ValueError("차수는 정수로 입력하세요.")
    if k <= 0:
        raise ValueError("차수는 양수여야 합니다.")

    matrix = []
    for i in range(k):
        while True:
            row_input = input(f"{i + 1}행 입력: ").strip()
            parts = row_input.split()
            if len(parts) != k:
                print(f"{k}개의 실수를 공백으로 구분해 입력하세요.")
                continue
            try:
                vals = [float(x) for x in parts]
                matrix.append(vals)
                break
            except ValueError:
                print("숫자만 입력하세요.")
    return matrix

# ---------------------------
# Main
# ---------------------------
def main():
    try:
        A = read_square_matrix()

        if not has_inverse(A):
            print("오류! 행렬식이 0이므로 역행렬을 구할 수 없습니다.")
            return

        print("\n[가우스-조던 소거법으로 구한 역행렬]\n")
        inv_gj = inverse_by_gauss_jordan(A)
        printMatrix(inv_gj)

        print("\n[행렬식(여인수/수반행렬)으로 구한 역행렬]\n")
        inv_det = inverse_by_determinant(A)
        printMatrix(inv_det)

        print("\n[결과 비교]")
        if are_same(inv_gj, inv_det, tol=1e-6):
            print("두 방법의 결과가 동일합니다. (tol=1e-6)")
        else:
            print("두 방법의 결과가 다릅니다. (수치 오차 또는 입력 확인)")

        print("\n[A * (가우스-조던 역행렬) ≈ I 확인]")
        prod1 = multiply(A, inv_gj)
        printMatrix(prod1)
        print("→ 단위행렬 여부:", "예" if is_identity(prod1) else "아니오")

        print("\n[A * (행렬식 역행렬) ≈ I 확인]")
        prod2 = multiply(A, inv_det)
        printMatrix(prod2)
        print("→ 단위행렬 여부:", "예" if is_identity(prod2) else "아니오")

    except ValueError as e:
        print(f"입력/계산 오류: {e}")

if __name__ == "__main__":
    main()
