# Jacobian-Curve-Libgcrypt

Jacobian-Curve-Libgcrypt - an implementation of the scalar multiplication method on the Jacobian elliptic curve using the doubling addition algorithm. 

## Author

It was created by Alexander Khamitov for the subject "Programming Algorithms of Protection of the Information" in HSE (2019-2020)

## Results

Test 1

Point E coordinates:
X: 00
Y: 01
Does the point lie on the curve? YES

Point E + E coordinates:
X: 00
Y: 01
Does the point lie on the curve? YES

Are points E and E + E equal? YES

-----------------------------------------------------------------------------

Test 2

Point B coordinates:
X: 1A
Y: 480C9BDE2C8F8F18BDF13E073A7776DD5ED2239AEFFB63116D84C2CA4DC6C432
Does the point lie on the curve? YES

Point E + B coordinates:
X: 1A
Y: 480C9BDE2C8F8F18BDF13E073A7776DD5ED2239AEFFB63116D84C2CA4DC6C432
Does the point lie on the curve? YES

Are points B and E + B equal? YES

Point B + E coordinates:
X: 1A
Y: 480C9BDE2C8F8F18BDF13E073A7776DD5ED2239AEFFB63116D84C2CA4DC6C432
Does the point lie on the curve? YES

Are points B and B + E equal? YES

-----------------------------------------------------------------------------

Test 3

Point B + B coordinates:
X: 00ED0E1ACF28915C7CA255D66B6001C8986B68EB1D39C14DDB66E9EFA9ADC127DB
Y: 19C18451736894BB8061D54C160D8E9588A72AE8F9BCCD42553A9B5BAB2F2C98
Does the point lie on the curve? YES

Point 2*B coordinates:
X: 00ED0E1ACF28915C7CA255D66B6001C8986B68EB1D39C14DDB66E9EFA9ADC127DB
Y: 19C18451736894BB8061D54C160D8E9588A72AE8F9BCCD42553A9B5BAB2F2C98
Does the point lie on the curve? YES

Are points B + B and 2*B equal? YES

-----------------------------------------------------------------------------

Test 4

Point (q + 1)*B coordinates:
X: 1A
Y: 480C9BDE2C8F8F18BDF13E073A7776DD5ED2239AEFFB63116D84C2CA4DC6C432
Does the point lie on the curve? YES

Are points (q + 1)*B and B equal? YES

-----------------------------------------------------------------------------

Test 5

Point q*B coordinates:
X: 00
Y: 01
Does the point lie on the curve? YES

Are points q*B and E equal? YES

-----------------------------------------------------------------------------

Test 6

Point (q - 1)*B coordinates:
X: 00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD7D
Y: 480C9BDE2C8F8F18BDF13E073A7776DD5ED2239AEFFB63116D84C2CA4DC6C432
Does the point lie on the curve? YES

Point -B coordinates:
X: -1A
Y: 480C9BDE2C8F8F18BDF13E073A7776DD5ED2239AEFFB63116D84C2CA4DC6C432
Does the point lie on the curve? YES

Are points (q - 1)*B and -B equal? YES

-----------------------------------------------------------------------------

Test 7

Point B + B + ... + B (scalar times) coordinates:
X: 00CC9BAD29888C0F06ABDB6074FB2C27684D100858EBFA4B899BD14C5E22ACDC77
Y: 00C7A9CCD9AA0240CE2D5E43BF08BB6B8F6354B86C6C8603DFCA66A51124DB0EFF
Does the point lie on the curve? YES

Point scalar*B coordinates:
X: 00CC9BAD29888C0F06ABDB6074FB2C27684D100858EBFA4B899BD14C5E22ACDC77
Y: 00C7A9CCD9AA0240CE2D5E43BF08BB6B8F6354B86C6C8603DFCA66A51124DB0EFF
Does the point lie on the curve? YES

Are points B + B + ... + B (957 times) and 957*B equal? YES

-----------------------------------------------------------------------------

Test 8

Point k1*B + k2*B coordinates:
X: 00BE15DFB7F61A2EFA8F3C5F8BEA4170D14AB24189CBE932561BF722C979C39E1F
Y: 009D78D46FA99C4DE24C77F345121B640E2B075DB153EA6139A461A5B443F1DFEE
Does the point lie on the curve? YES

Point (k1 + k2)*B:
X: 00BE15DFB7F61A2EFA8F3C5F8BEA4170D14AB24189CBE932561BF722C979C39E1F
Y: 009D78D46FA99C4DE24C77F345121B640E2B075DB153EA6139A461A5B443F1DFEE
Does the point lie on the curve? YES

Are points k1*B + k2*B and (k1 + k2)*B equal? YES

-----------------------------------------------------------------------------