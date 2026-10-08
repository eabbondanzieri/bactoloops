"""ImageJ's Process > Binary > Watershed, ported line for line from ImageJ 1.54p
(ij/plugin/filter/EDM.java: makeFloatEDM; ij/plugin/filter/MaximumFinder.java: findMaxima with
tolerance 0.5, output SEGMENTED, isEDM=true, strict=false, no edge exclusion).

watershed(binary) takes a boolean image (True = particle) and returns the segmented boolean image
ImageJ produces. Compiled with numba; the first call takes a few seconds to compile.
"""
import numpy as np
from numba import njit

NO_POINT = -1
MAXIMUM, LISTED, PROCESSED, MAX_AREA, EQUAL, MAX_POINT, ELIMINATED = 1, 2, 4, 8, 16, 32, 64
SQRT2 = np.float32(1.4142135624)
DX = np.array([0, 1, 1, 1, 0, -1, -1, -1], np.int64)
DY = np.array([-1, -1, 0, 1, 1, 1, 0, -1], np.int64)


# ------------------------------------------------------------------ EDM.makeFloatEDM
@njit(cache=True)
def _min_dist2(points, p_prev, p_diag, x, y, dist_sqr):
    p0 = points[x]
    nearest = p0
    if p0 != NO_POINT:
        x0 = p0 & 0xffff
        y0 = (p0 >> 16) & 0xffff
        d1 = (x - x0) * (x - x0) + (y - y0) * (y - y0)
        if d1 < dist_sqr:
            dist_sqr = d1
    if p_diag != p0 and p_diag != NO_POINT:
        x1 = p_diag & 0xffff
        y1 = (p_diag >> 16) & 0xffff
        d1 = (x - x1) * (x - x1) + (y - y1) * (y - y1)
        if d1 < dist_sqr:
            nearest = p_diag
            dist_sqr = d1
    if p_prev != p_diag and p_prev != NO_POINT:
        x1 = p_prev & 0xffff
        y1 = (p_prev >> 16) & 0xffff
        d1 = (x - x1) * (x - x1) + (y - y1) * (y - y1)
        if d1 < dist_sqr:
            nearest = p_prev
            dist_sqr = d1
    points[x] = nearest
    return np.float32(dist_sqr)


@njit(cache=True)
def _edm_line(fg, f, buf0, buf1, width, y):
    big = np.int64(2147483647)                    # Integer.MAX_VALUE: edges are not background
    p_prev = NO_POINT
    p_diag = NO_POINT
    for x in range(width):
        p_next_diag = buf0[x]
        if not fg[y, x]:
            buf0[x] = x | (y << 16)
        else:
            d2 = _min_dist2(buf0, p_prev, p_diag, x, y, big)
            if f[y, x] > d2:
                f[y, x] = d2
        p_prev = buf0[x]
        p_diag = p_next_diag
    p_prev = NO_POINT
    p_diag = NO_POINT
    for x in range(width - 1, -1, -1):
        p_next_diag = buf1[x]
        if not fg[y, x]:
            buf1[x] = x | (y << 16)
        else:
            d2 = _min_dist2(buf1, p_prev, p_diag, x, y, big)
            if f[y, x] > d2:
                f[y, x] = d2
        p_prev = buf1[x]
        p_diag = p_next_diag


@njit(cache=True)
def float_edm(fg):
    h, w = fg.shape
    f = np.zeros((h, w), np.float32)
    fmax = np.float32(3.4028234663852886e38)
    for y in range(h):
        for x in range(w):
            if fg[y, x]:
                f[y, x] = fmax
    buf0 = np.full(w, NO_POINT, np.int64)
    buf1 = np.full(w, NO_POINT, np.int64)
    for y in range(h):
        _edm_line(fg, f, buf0, buf1, w, y)
    buf0[:] = NO_POINT
    buf1[:] = NO_POINT
    for y in range(h - 1, -1, -1):
        _edm_line(fg, f, buf0, buf1, w, y)
    for y in range(h):
        for x in range(w):
            f[y, x] = np.float32(np.sqrt(f[y, x]))
    return f


# ------------------------------------------------------------------ MaximumFinder
@njit(cache=True)
def _within(x, y, d, w, h):
    xm, ym = w - 1, h - 1
    if d == 0:
        return y > 0
    if d == 1:
        return x < xm and y > 0
    if d == 2:
        return x < xm
    if d == 3:
        return x < xm and y < ym
    if d == 4:
        return y < ym
    if d == 5:
        return x > 0 and y < ym
    if d == 6:
        return x > 0
    return x > 0 and y > 0


@njit(cache=True)
def _true_edm_height(x, y, p, w, h, diro):
    off = x + y * w
    v = p[off]
    if x == 0 or y == 0 or x == w - 1 or y == h - 1 or v == 0:
        return v
    true_h = np.float32(v + np.float32(0.5) * SQRT2)
    ridge_or_max = False
    for d in range(4):
        d2 = (d + 4) % 8
        v1 = p[off + diro[d]]
        v2 = p[off + diro[d2]]
        if v >= v1 and v >= v2:
            ridge_or_max = True
            hh = np.float32((v1 + v2) / np.float32(2))
        else:
            hh = min(v1, v2)
        hh = np.float32(hh + (np.float32(1) if d % 2 == 0 else SQRT2))
        if true_h > hh:
            true_h = hh
    if not ridge_or_max:
        true_h = v
    return true_h


@njit(cache=True)
def _fate_table():
    table = np.zeros(256, np.int64)
    for item in range(256):
        is_set = np.zeros(8, np.bool_)
        for i in range(8):
            is_set[i] = (item >> i) & 1 == 1
        for i in range(8):
            if is_set[(i + 4) % 8]:
                table[item] |= 1 << i
        for i in range(0, 8, 2):
            if is_set[i]:
                is_set[(i + 1) % 8] = True
                is_set[(i + 7) % 8] = True
        transitions = 0
        for i in range(8):
            if is_set[i] != is_set[(i + 1) % 8]:
                transitions += 1
        if transitions >= 4:
            table[item] = 0
    return table


@njit(cache=True)
def _segment(p, w, h, tolerance):
    """findMaxima(ip=EDM, tolerance, NO_THRESHOLD, SEGMENTED, excludeOnEdges=false, isEDM=true)"""
    n = w * h
    diro = np.array([-w, -w + 1, 1, w + 1, w, w - 1, -1, -w - 1], np.int64)
    types = np.zeros(n, np.int64)
    gmin = np.float32(3.4028234663852886e38)
    gmax = np.float32(-3.4028234663852886e38)
    for i in range(n):
        if gmin > p[i]:
            gmin = p[i]
        if gmax < p[i]:
            gmax = p[i]
    out = np.zeros(n, np.int64)
    if not gmax > gmin:
        return out
    # getSortedMaxPoints
    n_max = 0
    for y in range(h):
        for x in range(w):
            i = x + y * w
            v = p[i]
            if v == gmin:
                continue
            v_true = _true_edm_height(x, y, p, w, h, diro)
            is_max = True
            inner = y != 0 and y != h - 1 and x != 0 and x != w - 1
            for d in range(8):
                if inner or _within(x, y, d, w, h):
                    vn = p[i + diro[d]]
                    vn_true = _true_edm_height(x + DX[d], y + DY[d], p, w, h, diro)
                    if vn > v and vn_true > v_true:
                        is_max = False
                        break
            if is_max:
                types[i] = MAXIMUM
                n_max += 1
    vfactor = np.float32(2e9 / (gmax - gmin))
    keys = np.empty(n_max, np.int64)
    k = 0
    for y in range(h):
        for x in range(w):
            i = x + y * w
            if types[i] == MAXIMUM:
                fv = _true_edm_height(x, y, p, w, h, diro)
                iv = np.int64(np.int32(np.float32(np.float32(fv - gmin) * vfactor)))
                keys[k] = (iv << 32) | i
                k += 1
    keys = np.sort(keys)
    # analyzeAndMarkMaxima (strict=false, excludeEdgesNow=false)
    max_sort_err = np.float32(np.float32(1.1) * (SQRT2 / np.float32(2)))
    tol = np.float32(tolerance)
    plist = np.empty(n, np.int64)
    for im in range(n_max - 1, -1, -1):
        off0 = keys[im] & 0xffffffff
        if types[off0] & PROCESSED:
            continue
        x0 = off0 % w
        y0 = off0 // w
        v0 = _true_edm_height(x0, y0, p, w, h, diro)
        sorting_error = True
        while sorting_error:
            plist[0] = off0
            types[off0] |= EQUAL | LISTED
            list_len = 1
            li = 0
            sorting_error = False
            max_possible = True
            x_eq = float(x0)
            y_eq = float(y0)
            n_eq = 1
            while True:
                off = plist[li]
                x = off % w
                y = off // w
                inner = y != 0 and y != h - 1 and x != 0 and x != w - 1
                for d in range(8):
                    off2 = off + diro[d]
                    if (inner or _within(x, y, d, w, h)) and (types[off2] & LISTED) == 0:
                        if p[off2] <= 0:
                            continue
                        if types[off2] & PROCESSED:
                            max_possible = False
                            break
                        x2 = x + DX[d]
                        y2 = y + DY[d]
                        v2 = _true_edm_height(x2, y2, p, w, h, diro)
                        if v2 > np.float32(v0 + max_sort_err):
                            max_possible = False
                            break
                        elif v2 >= np.float32(v0 - tol):
                            if v2 > v0:
                                sorting_error = True
                                off0 = off2
                                v0 = v2
                                x0 = x2
                                y0 = y2
                            plist[list_len] = off2
                            list_len += 1
                            types[off2] |= LISTED
                            if v2 == v0:
                                types[off2] |= EQUAL
                                x_eq += x2
                                y_eq += y2
                                n_eq += 1
                if not max_possible:
                    pass
                li += 1
                if li >= list_len:
                    break
                if not max_possible:
                    # Java breaks out of the direction loop only; the do-while continues
                    pass
            if sorting_error:
                for q in range(list_len):
                    types[plist[q]] = 0
            else:
                reset_mask = ~(LISTED if max_possible else (LISTED | EQUAL))
                x_eq /= n_eq
                y_eq /= n_eq
                min_d2 = 1e20
                nearest = 0
                for q in range(list_len):
                    off = plist[q]
                    x = off % w
                    y = off // w
                    types[off] &= reset_mask
                    types[off] |= PROCESSED
                    if max_possible:
                        types[off] |= MAX_AREA
                        if types[off] & EQUAL:
                            d2 = (x_eq - x) * (x_eq - x) + (y_eq - y) * (y_eq - y)
                            if d2 < min_d2:
                                min_d2 = d2
                                nearest = q
                if max_possible:
                    types[plist[nearest]] |= MAX_POINT
    if n_max == 0:
        for i in range(n):
            types[i] = PROCESSED | MAX_AREA
    # make8bit (isEDM: threshold 0.5, minValue 1)
    min_value = 1.0
    offset = min_value - (float(gmax) - min_value) * (1.0 / 253 / 2 - 1e-6)
    factor = 253.0 / (float(gmax) - min_value)
    if factor > 1:
        factor = 1.0
    pix = np.zeros(n, np.int64)
    for i in range(n):
        raw = p[i]
        if raw < 0.5:
            pix[i] = 0
        elif types[i] & MAX_AREA:
            pix[i] = 255
        else:
            v = 1 + np.int64(np.floor((float(raw) - offset) * factor + 0.5))
            if v < 1:
                pix[i] = 1
            elif v <= 254:
                pix[i] = v
            else:
                pix[i] = 254
    # cleanupMaxima
    for im in range(n_max - 1, -1, -1):
        off0 = keys[im] & 0xffffffff
        if types[off0] & (MAX_AREA | ELIMINATED):
            continue
        level = pix[off0]
        lo = level + 1
        plist[0] = off0
        types[off0] |= LISTED
        list_len = 1
        last_len = 1
        saddle = False
        while (not saddle) and lo > 0:
            lo -= 1
            last_len = list_len
            li = 0
            while True:
                off = plist[li]
                x = off % w
                y = off // w
                inner = y != 0 and y != h - 1 and x != 0 and x != w - 1
                for d in range(8):
                    off2 = off + diro[d]
                    if (inner or _within(x, y, d, w, h)) and (types[off2] & LISTED) == 0:
                        if (types[off2] & MAX_AREA) or ((types[off2] & ELIMINATED) and pix[off2] >= lo):
                            saddle = True
                            break
                        elif pix[off2] >= lo and (types[off2] & ELIMINATED) == 0:
                            plist[list_len] = off2
                            list_len += 1
                            types[off2] |= LISTED
                if saddle:
                    break
                li += 1
                if li >= list_len:
                    break
        for q in range(list_len):
            types[plist[q]] &= ~LISTED
        for q in range(last_len):
            pix[plist[q]] = lo
            types[plist[q]] |= ELIMINATED
    # watershedSegment
    hist = np.zeros(256, np.int64)
    for i in range(n):
        hist[pix[i]] += 1
    array_size = n - hist[0] - hist[255]
    coords = np.empty(max(array_size, 1), np.int64)
    highest = 0
    max_bin = 0
    level_start = np.zeros(256, np.int64)
    o = 0
    for v in range(1, 255):
        level_start[v] = o
        o += hist[v]
        if hist[v] > 0:
            highest = v
        if hist[v] > max_bin:
            max_bin = hist[v]
    level_off = np.zeros(highest + 1, np.int64)
    for i in range(n):
        v = pix[i]
        if v > 0 and v < 255:
            coords[level_start[v] + level_off[v]] = i
            level_off[v] += 1
    set_list = np.empty(max(min(max_bin, (n + 2) // 3), 1), np.int64)
    table = _fate_table()
    dir_seq = np.array([7, 3, 1, 5, 0, 4, 2, 6], np.int64)
    xmax, ymax = w - 1, h - 1
    for level in range(highest, 0, -1):
        remaining = hist[level]
        idle = 0
        while remaining > 0 and idle < 8:
            d_index = 0
            while True:
                pas = dir_seq[d_index % 8]
                # processLevel
                n_changed = 0
                n_unchanged = 0
                ls = level_start[level]
                for q in range(remaining):
                    off = coords[ls + q]
                    x = off % w
                    y = off // w
                    index = 0
                    if y > 0 and pix[off - w] == 255:
                        index ^= 1
                    if x < xmax and y > 0 and pix[off - w + 1] == 255:
                        index ^= 2
                    if x < xmax and pix[off + 1] == 255:
                        index ^= 4
                    if x < xmax and y < ymax and pix[off + w + 1] == 255:
                        index ^= 8
                    if y < ymax and pix[off + w] == 255:
                        index ^= 16
                    if x > 0 and y < ymax and pix[off + w - 1] == 255:
                        index ^= 32
                    if x > 0 and pix[off - 1] == 255:
                        index ^= 64
                    if x > 0 and y > 0 and pix[off - w - 1] == 255:
                        index ^= 128
                    mask = 1 << pas
                    if (table[index] & mask) == mask:
                        set_list[n_changed] = off
                        n_changed += 1
                    else:
                        coords[ls + n_unchanged] = off
                        n_unchanged += 1
                for q in range(n_changed):
                    pix[set_list[q]] = 255
                remaining -= n_changed
                if n_changed > 0:
                    idle = 0
                d_index += 1
                if not (remaining > 0):
                    break
                cond = idle < 8
                idle += 1
                if not cond:
                    break
        if remaining > 0 and level > 1:
            nxt = level
            while True:
                nxt -= 1
                if not (nxt > 1 and hist[nxt] == 0):
                    break
            if nxt > 0:
                new_end = level_start[nxt] + hist[nxt]
                pstart = level_start[level]
                for q in range(remaining):
                    off = coords[pstart + q]
                    x = off % w
                    y = off // w
                    add = False
                    if x == 0 or y == 0 or x == w - 1 or y == h - 1:
                        add = True
                    else:
                        for d in range(8):
                            if _within(x, y, d, w, h) and pix[off + diro[d]] == 0:
                                add = True
                                break
                    if add:
                        coords[new_end] = off
                        new_end += 1
                hist[nxt] = new_end - level_start[nxt]
    # watershedPostProcess
    for i in range(n):
        out[i] = 255 if pix[i] == 255 else 0
    return out


def watershed(binary):
    """Process > Binary > Watershed on a boolean image (True = particle)."""
    fg = np.ascontiguousarray(binary.astype(np.bool_))
    edm = float_edm(fg)
    h, w = fg.shape
    seg = _segment(edm.ravel(), w, h, 0.5)
    return seg.reshape(h, w) == 255
